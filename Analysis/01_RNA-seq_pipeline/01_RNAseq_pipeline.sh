#!/usr/bin/bash

help() {
		echo ""
		echo -e "Usage: \n\t bash $0 [options] -d <project_dir> -o <output_dir> -i <data_dir> --ref <reference_genome> --gtf <GTF_flie>"
		echo ""
		echo -e "  -d \n\t\t Project directory."
		echo -e "  -i \n\t\t Data directory for raw and trimmed fastqs. All raw fastqs should be put in ${data_dir}/fastq."
		echo -e "  -o \n\t\t Output directory for results."
		echo -e "  --ref \n\t\t Reference for STAR mapping."
		echo -e "  --gtf \n\t\t GTF for gene quantification."
		echo ""
		echo "Options:"
		echo -e "  -t, --thread \n\t\t Numbers of threads, default: 1"
		echo -e "  --pair \n\t\t Flag to turn on pair-end mode."
		echo -e "  -h, --help \n\t\t show help message."
		echo ""
}

if [ $# -eq 0 ]; then
        echo "Run 'bash $0 -h' to see more information."
        exit 0
fi

# NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
# separately; see below.
TEMP=`getopt -o d:o:t:i:h --long help,ref:,threads:,pair,gtf: \
             -n 'RNAseq_pipeline.sh' -- "$@"`

if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi

# Note the quotes around `$TEMP': they are essential!
eval set -- "$TEMP"

THREADS=1
PAIR=false
FQ1=
FQ2=
# Illumina universial adapters
#FORADAPT=AGATCGGAAGAG
#READAPT=AGATCGGAAGAG


while true; do
  case "$1" in
				-d ) PROJECTDIR="$2"; shift 2;;
				-o ) OUTDIR="$2"; shift 2;;
				-i ) DATADIR="$2"; shift 2;;
				-t | --threads ) THREADS="$2"; shift 2;;
				--pair ) PAIR=true; shift ;;
				--ref ) REFERENCE="$2"; shift 2;;
				--gtf ) GTF="$2"; shift 2;;
				-h | --help)
                        help
                        exit 0;;
                -- ) shift; break ;;
                * ) echo "Invalid option: ${optionName}" ;
				   echo "Try 'bash `basename $0` -h' for more information" ; 
				   break ;;
        esac
done

if [[ ! -d ${PROJECTDIR}/data ]] || [[ ! -d ${PROJECTDIR}/results ]]; then
	echo " The 'data' and 'results' folders must be in your project directory"
	exit 1
fi

if [[ ! -d ${OUTDIR} ]]; then
	echo " The output directory needed to be specified"
	exit 1
fi

if [[ ! -d ${DATADIR} ]]; then
	echo " The data directory needed to be specified"
	exit 1
fi

if [[ $REFERENCE = "" ]]; then
	echo " The reference genome is needed "
	exit 1
fi

## mkdir
if [[ ! -d ${DATADIR}/clean ]]; then mkdir -p ${DATADIR}/clean; fi
FASTQLOC=${DATADIR}/fastq
CLEANLOC=${DATADIR}/clean
if [[ ! -d ${OUTDIR}/align ]]; then mkdir -p ${OUTDIR}/align ; fi
ALIGNLOC=${OUTDIR}/align
if [[ ! -d ${OUTDIR}/featurecounts ]]; then mkdir -p ${OUTDIR}/featurecounts ; fi
COUNTLOC=${OUTDIR}/featurecounts
if [[ ! -d ${OUTDIR}/QC ]]; then mkdir -p ${OUTDIR}/QC; fi

## Printing 
echo ""
echo " Project directory: ${PROJECTDIR} "
echo " Fastq directory: ${FASTQLOC}"
echo " Reference path: ${REFERENCE}"
echo " Numbers of threads: ${THREADS} "
echo ""
## list all fastq files
FQs=(`ls ${FASTQLOC}/*fastq.gz`)
## Processing with for loop
if $PAIR; then
	## Paired-end mode
	for ((j=0; j<=${#FQs[@]}-1; j+=2)); do
		FQ1=`basename ${FQs[$j]}`
		FQ2=`basename ${FQs[$j+1]}`
		echo "Processing " ${FQ1}" "${FQ2}
		## Extract the common prefix of fastq1 and fastq2
		PREFIX=`{ echo "$FQ1"; echo "$FQ2";  } | sed -e 'N;s/^\(.*\).*\n\1.*$/\1\n\1/;D'`

		## trimming
		echo ""
		echo ">Step01 Trimming & QC"
		echo ""
		trim_galore --nextseq 30 \
		--phred33 \
		--gzip \
		-o ${CLEANLOC}/ \
		--cores ${THREADS} \
		--basename ${PREFIX} \
		--fastqc_args "-o ${OUTDIR}/QC -t 8" \
		--paired \
		${FASTQLOC}/${FQ1} ${FASTQLOC}/${FQ2}

		## Mapping -- ../align/{PREFIX}/
		echo ""
		echo ">Step02 Mapping" 
		echo ""
		if [[ ! -d ${ALIGNLOC}/${PREFIX}_align ]]; then mkdir -p ${ALIGNLOC}/${PREFIX}_align; fi
		STAR --genomeDir ${REFERENCE} \
		--readFilesCommand zcat \
		--readFilesIn ${CLEANLOC}/${PREFIX}_val_1.fq.gz ${CLEANLOC}/${PREFIX}_val_2.fq.gz \
		--runThreadN ${THREADS} \
		--outSAMtype SAM \
		--outFileNamePrefix ${ALIGNLOC}/${PREFIX}_align/ 

		## samtools -- ../align/{PREFIX}/
		echo ""
		echo ">Step03 SAMTOOLS" 
		echo ""
		samtools view -@ ${THREADS} -q 30 -f 2 -hSb ${ALIGNLOC}/${PREFIX}_align/Aligned.out.sam |samtools sort - -@ ${THREADS} -o ${ALIGNLOC}/${PREFIX}_align/${PREFIX}.sorted.bam
		samtools index ${ALIGNLOC}/${PREFIX}_align/${PREFIX}.sorted.bam
		rm -f ${ALIGNLOC}/${PREFIX}_align/Aligned.out.sam
		
		## featureCounts -- ${OUTDIR}/featurecounts
		echo ""
		echo ">Step04 featureCounts" 
		echo ""
		featureCounts \
		-p -B -C \
		-t exon \
		-g gene_id \
		-T ${THREADS} \
		-a $GTF \
		-o ${COUNTLOC}/${PREFIX}_counts.txt \
		${ALIGNLOC}/${PREFIX}_align/${PREFIX}.sorted.bam
	done
else
	## Single-end mode
	for ((j=0; j<=${#FQs[@]}-1; j++)); do
		FQ=`basename ${FQs[$j]}`
		echo "Processing " ${FQ}
		PREFIX=`basename ${FQ} .fastq.gz`

		## trimming
		echo ""
		echo ">Step01 Trimming & QC"
		echo ""
		trim_galore --nextseq 30 \
		--phred33 \
		--gzip \
		-o ${CLEANLOC}/ \
		--cores ${THREADS} \
		--basename ${PREFIX} \
		--fastqc_args "-o ${OUTDIR}/QC -t 8" \
		${FASTQLOC}/${FQ} 

		## Mapping -- ../align/{PREFIX}/
		echo ""
		echo ">Step02 Mapping" 
		echo ""
		if [[ ! -d ${ALIGNLOC}/${PREFIX}_align ]]; then mkdir -p ${ALIGNLOC}/${PREFIX}_align; fi
		STAR --genomeDir ${REFERENCE} \
		--readFilesCommand zcat \
		--readFilesIn ${CLEANLOC}/${PREFIX}_trimmed.fq.gz \
		--runThreadN ${THREADS} \
		--outSAMtype SAM \
		--outFileNamePrefix ${ALIGNLOC}/${PREFIX}_align/  

		## samtools -- ../align/{PREFIX}/
		echo ""
		echo ">Step03 SAMTOOLS" 
		echo ""
		samtools view -@ ${THREADS} -q 30 -F 4 -hSb ${ALIGNLOC}/${PREFIX}_align/Aligned.out.sam |samtools sort - -@ ${THREADS} -o ${ALIGNLOC}/${PREFIX}_align/${PREFIX}.sorted.bam
		samtools index ${ALIGNLOC}/${PREFIX}_align/${PREFIX}.sorted.bam
		rm -f ${ALIGNLOC}/${PREFIX}_align/Aligned.out.sam

		## featureCounts -- ${OUTDIR}/featurecounts
		echo ""
		echo ">Step04 featureCounts" 
		echo ""
		featureCounts \
		-t exon \
		-g gene_id \
		-T ${THREADS} \
		-a $GTF \
		-o ${COUNTLOC}/${PREFIX}_counts.txt \
		${ALIGNLOC}/${PREFIX}_align/${PREFIX}.sorted.bam
	done
fi


# post multiqc
## trimmed reads qc 
multiqc -o ${OUTDIR}/QC/  ${OUTDIR}/QC/*zip
## star alignment reports
multiqc -o ${OUTDIR}/align/  ${OUTDIR}/align/*/Log.final.out
