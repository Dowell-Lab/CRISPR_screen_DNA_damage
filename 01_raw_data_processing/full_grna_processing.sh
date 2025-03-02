#!/bin/bash

seqrun=seq1
fastq=/path/to/fastq
script=/path/to/scripts
process=/path/to/general/processing/files
scratch=/path/to/scratch/"$seqrun"

mkdir -p "$scratch" "$scratch"/index
rsync "$process"/nt_guides_and_barcodes.txt "$scratch"/
rsync "$process"/grna_hisat/* "$scratch"/index/  # HISAT2 indeces for grna_library.fa

# Run sbatch script for each R1 fastq file
for sampfile in "$fastq"/*"$seqrun"*.fastq.gz; do

sbatch --export=sampfile=$sampfile,seqrun=$seqrun "$script"/full_grna_processing.sbatch

done
