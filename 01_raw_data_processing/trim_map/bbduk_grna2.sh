#!/bin/bash

# Inputs

scratch=$1
filename=$2
filebase="$1"/trim/"$3"
ncores=$4

# Process Read 2 for gRNA2

bbduk.sh -Xmx100g \
    t="$ncores" \
    in="$scratch"/trim/"$filename" \
    out="$filebase"_R2.unmatched_trim1.fastq.gz \
    literal=TTTCGCTGACGT \
    k=12 hdist=1 rcomp=f \
    stats="$filebase"_stats1.txt

bbduk.sh -Xmx100g \
    t="$ncores" \
    in="$filebase"_R2.unmatched_trim1.fastq.gz \
    out="$filebase"_R2.unmatched_trim2.fastq.gz \
    outm="$filebase"_R2.matched_trim2.fastq.gz \
    literal=ATTTCTAGCTCTAAAAC \
    k=17 hdist=2 rcomp=f \
    stats="$filebase"_stats2.txt
rm "$filebase"_R2.unmatched_trim1.fastq.gz
rm "$filebase"_R2.unmatched_trim2.fastq.gz

bbduk.sh -Xmx100g \
    t="$ncores" \
    in="$filebase"_R2.matched_trim2.fastq.gz \
    out="$filebase"_R2.unmatched_trim3.fastq.gz \
    outm="$filebase"_R2.matched_trim3.fastq.gz \
    literal=CGGTGTTTCGTC \
    k=12 hdist=1 rcomp=f \
    stats="$filebase"_stats3.txt
rm "$filebase"_R2.matched_trim2.fastq.gz
rm "$filebase"_R2.unmatched_trim3.fastq.gz

bbduk.sh -Xmx100g \
    t="$ncores" \
    in="$filebase"_R2.matched_trim3.fastq.gz \
    out="$filebase"_R2.trim4.fastq.gz \
    literal=ATTTCTAGCTCTAAAAC \
    ktrim=l k=17 mink=8 hdist=2 rcomp=f \
    stats="$filebase"_stats4.txt
rm "$filebase"_R2.matched_trim3.fastq.gz

bbduk.sh -Xmx100g \
    t="$ncores" \
    in="$filebase"_R2.trim4.fastq.gz \
    out="$filebase"_grna2.fastq.gz \
    literal=CGGTGTTTCGTC \
    ktrim=r k=12 hdist=1 maq=10 rcomp=f \
    minlen=20 maxlen=20 \
    stats="$filebase"_stats5.txt
rm "$filebase"_R2.trim4.fastq.gz

echo "Reverse primer filter" > "$filebase"_R2_trimstats.txt
cat "$filebase"_stats1.txt >> "$filebase"_R2_trimstats.txt
echo "Before guide filter" >> "$filebase"_R2_trimstats.txt
cat "$filebase"_stats2.txt >> "$filebase"_R2_trimstats.txt
echo "After guide filter" >> "$filebase"_R2_trimstats.txt
cat "$filebase"_stats3.txt >> "$filebase"_R2_trimstats.txt
echo "Before guide trim" >> "$filebase"_R2_trimstats.txt
cat "$filebase"_stats4.txt >> "$filebase"_R2_trimstats.txt
echo "After guide trim" >> "$filebase"_R2_trimstats.txt
cat "$filebase"_stats5.txt >> "$filebase"_R2_trimstats.txt

rm "$filebase"_stats*
