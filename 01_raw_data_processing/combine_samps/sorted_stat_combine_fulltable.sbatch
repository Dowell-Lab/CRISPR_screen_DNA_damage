#!/bin/bash
#SBATCH --job-name=stat_combine  # Job name
#SBATCH --mail-type=NONE # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<mail> # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=1 # Number of CPU (processer cores i.e. tasks)
#SBATCH --time=4:00:00 # Time limit hrs:min:sec
#SBATCH -p short
#SBATCH --mem=50gb # Memory limit
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

module load python/3.6.3

#################################################################

# Set script variables
script=/path/to/scripts
output=/path/to/accum/results
scratch=/path/to/scratch
indir="$scratch"/output/ # Reading straight in from sorted_stat_combine scripts
outfilebase="$scratch"/full_output/fullbasal_
#outfilebase="$scratch"/full_output/IR_ # IR processed into different countfile
suffix=grna1_coupled_counts.tsv

#################################################################
mkdir -p "$scratch"/full_output "$output"

# Create combined files
touch "$outfilebase".tsv

# Parse and combine files
python3 "$script"/combine_samps/sorted_stat_combine_fulltable.py \
  "$indir" \
  "$outfilebase" \
  "$suffix"
  
rsync "$outfilebase"* "$output"/
