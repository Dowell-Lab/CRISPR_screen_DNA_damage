#!/bin/bash
#SBATCH --job-name=hisat  # Job name
#SBATCH --mail-type=NONE # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<email> # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=8 # Number of CPU (processer cores i.e. tasks)
#SBATCH --time=1:00:00 # Time limit hrs:min:sec
#SBATCH -p short
#SBATCH --mem=10gb # Memory limit
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

module load hisat2/2.1.0

#################################################################

# Set script variables
outdir=/path/to/processing/files/


#################################################################

mkdir -p "$outdir"/hisat_indeces

# Index genome

hisat2-build "$outdir"/grna_library.fa "$outdir"/hisat_indeces/grna

