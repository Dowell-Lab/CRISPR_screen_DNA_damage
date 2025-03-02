#!/bin/bash

# Inputs

scratch=$1
filename=$2
filebase=$3
seqtype=$4
ncores=$5

# Map gRNA1

if [ "$seqtype" = "grna1" ]
then
  hisat2 -p "$ncores" \
    -q \
    --ignore-quals \
    --norc \
    --no-softclip \
    --mp 2,1 \
    --score-min L,0,-0.4 \
    --rdg 20,20 \
    --rfg 20,20 \
    --no-unal \
    --no-spliced-alignment \
    -x "$scratch"/index/grna \
    -U "$filename" \
    --new-summary \
    > "$scratch"/sams/"$filebase"_"$seqtype".sam \
    2> "$scratch"/mapstats/"$filebase"_"$seqtype".mapstats.txt
else
  hisat2 -p "$ncores" \
    -q \
    --ignore-quals \
    --nofw \
    --no-softclip \
    --mp 2,1 \
    --score-min L,0,-0.4 \
    --rdg 20,20 \
    --rfg 20,20 \
    --no-unal \
    --no-spliced-alignment \
    -x "$scratch"/index/grna \
    -U "$filename" \
    --new-summary \
    > "$scratch"/sams/"$filebase"_"$seqtype".sam \
    2> "$scratch"/mapstats/"$filebase"_"$seqtype".mapstats.txt  
fi

rm "$filename"

# Cut SAM file to what I really need
cut -f 1,3,12,13 "$scratch"/sams/"$filebase"_"$seqtype".sam \
  | grep -v @ \
  > "$scratch"/sams/"$filebase"_"$seqtype".readinfo

rm "$scratch"/sams/"$filebase"_"$seqtype".sam
