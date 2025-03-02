#!/bin/bash
#SBATCH --job-name=fastqc  # Job name
#SBATCH --mail-type=NONE # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<mail> # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=1 # Number of CPU (processer cores i.e. tasks)
#SBATCH --time=12:00:00 # Time limit hrs:min:sec
#SBATCH -p short
#SBATCH --mem=8gb # Memory limit
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

module load fastqc/0.11.8

#################################################################

# Set script variables 
indir=/path/to/fastq
outdir=/path/to/outdir

#################################################################

mkdir -p "$outdir"/"$fqdate"

# Run a job on all fastq files in the directory
fastqc "$indir"/"$fqdate"/*{.fastq,.fastq.gz} -o "$outdir"/"$fqdate"

