#!/bin/bash
# Set default values for options
in_fasta=""
out=""

# Parse command line options
while [[ $# -gt 0 ]]
do
    key="$1"

    case $key in
        -i)
            in_fasta="$2"
            shift # past argument
            shift # past value
            ;;
        -o)
            out="$2"
            shift # past argument
            shift # past value
            ;;
        *)
            echo "Usage: $0 -i <input fasta> -o <output csv>"
            exit 1
            ;;
    esac
done

# make tempoary directory 
if [ ! -d tmp/ ];then
    mkdir tmp/
fi

# Read file line by line and echo its content two lines at a time
while read -r line1 && read -r line2; do
    echo -e "$line1\n$line2" > seqi.fa
    prefix=`echo $line1 |sed 's/>//' -`
    RNAfold -d2 --noLP < seqi.fa > tmp/${prefix}.out
done < "$in_fasta"
rm -rf *ps

# extract MEF from tmp/*.out files
if [ ! -f $out ]; then 
    touch $out
    echo -e "ID,MEF,coord" > $out
fi

for i in tmp/*.out; do
    id=`grep -o "ENSG[0-9]*" $i`
    mfe=`sed -n '3s/.*(\(.*\)).*/\1/p' $i`
    coord=`grep -o "chr.*" $i`
    echo -e "$id,$mfe,$coord" >> $out
done
rm -rf tmp/