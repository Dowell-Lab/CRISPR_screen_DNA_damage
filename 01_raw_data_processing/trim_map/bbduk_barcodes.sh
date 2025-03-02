#!/bin/bash

# Inputs

scratch=$1
filename=$2
filebase="$1"/trim/"$3"
ncores=$4

# Process Read 2 for barcodes

bbduk.sh -Xmx100g \
    t="$ncores" \
    in="$scratch"/trim/"$filename" \
    out="$filebase"_R2.unmatched_trim1_bc.fastq.gz \
    outm="$filebase"_R2.matched_trim1_bc.fastq.gz \
    literal=GCCCTACTATCCAATTG rcomp=f \
    k=17 hdist=2 \
    stats="$filebase"_stats1.txt
rm "$scratch"/trim/"$filename"
rm "$filebase"_R2.unmatched_trim1_bc.fastq.gz

bbduk.sh -Xmx100g \
    t="$ncores" \
    in="$filebase"_R2.matched_trim1_bc.fastq.gz \
    out="$filebase"_R2.unmatched_trim2_bc.fastq.gz \
    outm="$filebase"_R2.matched_trim2_bc.fastq.gz \
    literal=GAATTCCGTTGC \
    k=12 hdist=1 rcomp=f \
    stats="$filebase"_stats2.txt
rm "$filebase"_R2.matched_trim1_bc.fastq.gz
rm "$filebase"_R2.unmatched_trim2_bc.fastq.gz

bbduk.sh -Xmx100g \
    t="$ncores" \
    in="$filebase"_R2.matched_trim2_bc.fastq.gz \
    out="$filebase"_R2.trim3_bc.fastq.gz \
    literal=GCCCTACTATCCAATTG \
    ktrim=l k=17 mink=8 hdist=2 rcomp=f \
    stats="$filebase"_stats3.txt
rm "$filebase"_R2.matched_trim2_bc.fastq.gz

bbduk.sh -Xmx100g \
    t="$ncores" \
    in="$filebase"_R2.trim3_bc.fastq.gz \
    out="$filebase"_R2.trim4_bc.fastq.gz \
    literal=GAATTCCGTTGC \
    ktrim=r k=12 hdist=1 maq=10 rcomp=f \
    minlen=46 maxlen=46 \
    stats="$filebase"_stats4.txt
rm "$filebase"_R2.trim3_bc.fastq.gz

bbduk.sh -Xmx100g \
    t="$ncores" \
    in="$filebase"_R2.trim4_bc.fastq.gz \
    out="$filebase"_barcode2.fastq.gz \
    literal=GAATTG \
    ktrim=l k=6 hdist=0 maq=10 rcomp=f\
    minlen=20 maxlen=20 \
    stats="$filebase"_stats5_2.txt

bbduk.sh -Xmx100g \
    t="$ncores" \
    in="$filebase"_R2.trim4_bc.fastq.gz \
    out="$filebase"_barcode1.fastq.gz \
    literal=GAATTG \
    ktrim=r k=6 hdist=0 maq=10 rcomp=f \
    minlen=20 maxlen=20 \
    stats="$filebase"_stats5_1.txt
rm "$filebase"_R2.trim4_bc.fastq.gz

echo "Before barcode filter" > "$filebase"_barcode_trimstats.txt
cat "$filebase"_stats1.txt >> "$filebase"_barcode_trimstats.txt
echo "After barcode filter" >> "$filebase"_barcode_trimstats.txt
cat "$filebase"_stats2.txt >> "$filebase"_barcode_trimstats.txt
echo "Before barcode trim" >> "$filebase"_barcode_trimstats.txt
cat "$filebase"_stats3.txt >> "$filebase"_barcode_trimstats.txt
echo "After barcode trim" >> "$filebase"_barcode_trimstats.txt
cat "$filebase"_stats4.txt >> "$filebase"_barcode_trimstats.txt
echo "After barcode split:1" >> "$filebase"_barcode_trimstats.txt
cat "$filebase"_stats5_1.txt >> "$filebase"_barcode_trimstats.txt
echo "After barcode split:2" >> "$filebase"_barcode_trimstats.txt
cat "$filebase"_stats5_2.txt >> "$filebase"_barcode_trimstats.txt
rm "$filebase"_stats*
