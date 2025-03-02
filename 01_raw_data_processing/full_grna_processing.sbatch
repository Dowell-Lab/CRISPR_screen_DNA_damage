#!/bin/bash
#SBATCH --job-name=grna_process  # Job name
#SBATCH --mail-type=NONE # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<email> # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=64 # Number of CPU (processer cores i.e. tasks)
#SBATCH --time=10:00:00 # Time limit hrs:min:sec
#SBATCH -p short
#SBATCH --mem=100gb # Memory limit
#SBATCH --output=/path/to/outerr/%x_%j.out
#SBATCH --error=/path/to/outerr/%x_%j.err

################## JOB INFO #####################################

printf "\nDirectory: $INDIR"
printf "\nRun on: $(hostname)"
printf "\nRun from: $(pwd)"
printf "\nScript: $0\n"
date

printf "\nYou've requested $SLURM_CPUS_ON_NODE core(s).\n"

################## LOAD MODULES #################################

module load bbmap/38.05
module load fastqc/0.11.8
module load hisat2/2.1.0
module load python/3.6.3

#################################################################

# Set script variables
# When run in loop with full_grna_processing.sh, $sampfile and $seqrun are set there
#seqrun=seq1
#sampfile=/path/to/fastq/day0_rep1_"$seqrun"_R1.fastq.gz
ncores=64
script=/path/to/scripts
fastq=/path/to/fastq
output=/path/to/output/"$seqrun"
scratch=/path/to/scratch/"$seqrun"

#################################################################
mkdir -p "$scratch" "$scratch"/sams "$scratch"/trim \
  "$scratch"/mapstats "$scratch"/results "$scratch"/index \
  "$output"/trimstats "$output"/mapstats "$output"/results

filebase=$(basename "$sampfile" _R1.fastq.gz)
filename1=$(basename "$sampfile")
filename2=$(echo "$filename1" | sed -r 's/R1/R2/g')

rsync "$fastq"/"$filebase"* "$scratch"/trim/

# Trimming

sh "$script"/trim_map/bbduk_grna1.sh "$scratch" "$filename1" "$filebase" "$ncores"
sh "$script"/trim_map/bbduk_grna2.sh "$scratch" "$filename2" "$filebase" "$ncores"
sh "$script"/trim_map/bbduk_barcodes.sh "$scratch" "$filename2" "$filebase" "$ncores"

rsync "$scratch"/trim/"$filebase"_*_trimstats.txt "$output"/trimstats/
rm "$scratch"/trim/"$filebase"_*_trimstats.txt

# HISAT2 alignment

sh "$script"/trim_map/hisat_allseqs.sh \
  "$scratch" \
  "$scratch"/trim/"$filebase"_grna1.fastq.gz \
  "$filebase" \
  grna1 \
  "$ncores"

sh "$script"/trim_map/hisat_allseqs.sh \
  "$scratch" \
  "$scratch"/trim/"$filebase"_grna2.fastq.gz \
  "$filebase" \
  grna2 \
  "$ncores"

sh "$script"/trim_map/hisat_allseqs.sh \
  "$scratch" \
  "$scratch"/trim/"$filebase"_barcode1.fastq.gz \
  "$filebase" \
  barcode1 \
  "$ncores"

sh "$script"/trim_map/hisat_allseqs.sh \
  "$scratch" \
  "$scratch"/trim/"$filebase"_barcode2.fastq.gz \
  "$filebase" \
  barcode2 \
  "$ncores"

# Do read sorting

python3 "$script"/read_sorting/sort_seqs.py \
  "$scratch"/sams/ \
  "$filebase" \
  "$scratch"/results/ \
  "$scratch"/nt_guides_and_barcodes.txt
rm "$scratch"/sams/"$filebase"*

rsync "$scratch"/mapstats/"$filebase"* "$output"/mapstats/
rsync "$scratch"/results/"$filebase"* "$output"/results/
rm "$scratch"/results/"$filebase"*
