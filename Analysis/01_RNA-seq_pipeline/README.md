# Description

The `RNA_pipeline.sh`  integrates the common steps for RNA-seq analysis, including : 

1. Trimming and quality control (`Trim Galore` and `FastQC`)

2. Mapping (`STAR`)

3. sam2bam, sort and index (`samtools`)

4. gene expression quantification (`featureCounts`)




# Dependences

```shell
python #v3.7
fastqc #v0.11.9
trim_galore #v0.6.6
STAR #2.7.6a
samtools #1.11
featureCounts #2.0.1
```



# Prerequisite 

In order to perform analysis with `RNAseq_pipeline.sh`, a project directory and the sub-directories must be created as the following structure: 

```shell
<project_dir>/
|-- data
|   -- fastq
|   	`*fastq.gz`
|-- results
|-- src
```

Using the following codes to create required directories: 

```bash
# In <project_dir>
mkdir -p data/fastq
mkdir results
mkdir src
```

- `data/fastq/` for holding all `fastq` files to be analysis. Put all your `fastq` files in `<project_dir>/<data_dir>/fastq/`. 

- `results/` for holding results.

- `src/` for any source files, e.g., `RNAseq_pipeline.sh`..



# Usage
```Shell
bash RNAseq_pipeline.sh [--pair] -d <project_dir> -o <output_dir> -i <input_data_dir> --ref <reference_genome> --gtf <GTF_file> -t <threads>
```

`-d`: project directory, e.g.,  `~/RNA-seq_pipe`, absolute path is preferable. 

`-i`: data directory, un-processed `fastq` files are in `<data_dir>/fastq` and trimmed `fastq` files will be generated at `<data_dir>/clean`.

`-o`: output directory, create results in the specified output directory. 

`--ref`: path to the directory where genome files are stored. 

`--gtf`:  path to the annotation files (in GTF/GFF format).  Gzipped file is also accepted.

`--pair`: `FLAG` turn on paired-end processing mode. 

`-t`: `INT`, number of threads



## In enONE-paper

First, create a project directory, for example `RNApipeTest` 

```bash
mkdir RNApipe
```

In `<project_dir>`, create the directories mentioned above and copy all the `fastq` files in `<data_dir>/fastq` . Also, place `STAR_pipeline.sh` at `src/`. 

```bash
cd RNApipe
mkdir -p data/fastq
mkdir results
mkdir src
cp /wherever/you/store/the/fastq.gz ./data/fastq
```

After preparation of the project directory, it should look like:

```shell
RNApipe/
├── data
│   └── fastq
│       ├── S1_R1.fastq.gz
│       ├── S1_R2.fastq.gz
│       ├── S2_R1.fastq.gz
│       ├── S2_R2.fastq.gz
│       ├── ...
│       └── S20_R2.fastq.gz
├── results
└── src
 └── RNAseq_pipeline.sh
```

----



In enONE manuscript, we used the following commands. 

We first generated a combined genome and associated STAR index as well as a combined annotation file of human and fly. 

```bash
# pseudo codes
# combine annotation
cat Homo_sapiens.GRCh38.94.chr.gff3 dmel-all-r6.36.gff3 > GRCh38_Dm6.36.gff3

# combine genome
cat GRCh38.primary_assembly.genome.fa dmel-all-chromosome-r6.36.fasta > GRCh38_Dm6.36.fa

# generate index for combined genome
STAR --runMode genomeGenerate --genomeFastaFiles GRCh38_Dm6.36.fa --sjdbGTFfile GRCh38_Dm6.36.gff3 --sjdbOverhang 149 --runThreadN 16 --genomeDir GRCh38_Dm6.36/ --genomeSAindexNbases 12
```



> Genome:
>
> - human: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh38.primary_assembly.genome.fa.gz
> - fly: http://ftp.flybase.net/releases/FB2020_05/dmel_r6.36/fasta/dmel-all-chromosome-r6.36.fasta.gz
>
> Annotation:
>
> - human: https://ftp.ensembl.org/pub/release-94/gff3/homo_sapiens/Homo_sapiens.GRCh38.94.gff3.gz
> - fly: http://ftp.flybase.net/releases/FB2020_05/dmel_r6.36/gff/dmel-all-r6.36.gff.gz

Then, we used the combined references for sequence alignment and gene quantification. 

```shell
# In <project_dir>
cd src/
nohup bash RNAseq_pipeline.sh --pair \
-d ~/RNApipe \
-i ~/RNApipe/data \
-o ~/RNApipe/results \
--ref /Reference/hg_dm/index/star_index/GRCh38_Dm6.36/ \
--gtf /Reference/hg_dm/annotation/GRCh38_Dm6.36.gff3 \
-t 32 >nohup1.out 2>&1 &
```

Messages from the program will be directed to file `nohup1.out`, you can use `cat nohup1.out` or `tail nohup1.out` to check. 

 

# Output

1. The trimmed `fastqs` will be output in `<data_dir>/clean` 

```shell
data/clean/
|-- *_trimmed.fq.gz
...
```

2. The QC results of trimmed reads are in ` <output_dir>/QC/`

```shell
results/QC/
|-- *_R1_clean_fastqc.html
|-- *_R1_clean_fastqc.zip

```

Also, `multiqc` reports of `FastQC` will be generated at `results/QC/`. 

3.  The mapping results in `<output_dir>/star/`

```shell
results/star/
|-- *_align
|   |-- Aligned.out.sam
|   |-- Log.final.out
|   |-- Log.out
|   |-- Log.progress.out
|   |-- *.sorted.bam
|   |-- *.sorted.bam.bai
|   -- SJ.out.tab
...

```

Also, `multiqc` reports of `Log.final.out` from `STAR` alignment will be generated at`results/star/`. 

4. The gene expression quantification results in `<output_dir>/featurecounts/`

```shell
  results/featurecounts/
|-- *_counts
|   |-- *_counts.txt
|   `-- *_counts.txt.summary
```

Counts files in `results/featurecounts/*_counts.txt` can be used for further analysis. 



# Citations

If you use `RNAseq_pipeline.sh` for your analysis, please cite the publication: [Epitranscriptome analysis of NAD-capped RNA by spike-in-based normalization and prediction of chronological age]()




- **Author**: Dean Li
- **Date**: 2022-12-23

