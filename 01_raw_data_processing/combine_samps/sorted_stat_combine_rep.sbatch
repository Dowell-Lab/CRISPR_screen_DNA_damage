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
# When run in loop with sorted_stat_combine.sh, $sampid and $baseout are set there
#sampid=day0_rep1
#baseout=/path/to/read/count/output/parent/dir
script=/path/to/scripts
output=/path/to/final/accum/results
scratch=/path/to/scratch
infilebase="$scratch"/orig/"$sampid"
gtype=grna1
outfilebase="$scratch"/output/"$sampid"_"$gtype"

#################################################################
mkdir -p "$scratch" "$scratch"/orig "$scratch"/output "$output"

# Create combined files
touch "$outfilebase"_coupled.tsv
touch "$outfilebase"_fulldata.tsv

# Parse and combine files
for countfile in "$baseout"/*/results/"$sampid"*"$gtype"_coupled.tsv; do
  rsync "$countfile" "$infilebase"_orig_file.tsv
  awk '{print $2"\t"$3}' "$infilebase"_orig_file.tsv \
    | tail -n +2 \
    >> "$infilebase"_coupled.tsv
  rm "$infilebase"_orig_file.tsv
done

for countfile in "$baseout"/*/results/"$sampid"*"$gtype"_fulldata.tsv; do
  rsync "$countfile" "$infilebase"_orig_file.tsv
  awk '{print $2"\t"$3}' "$infilebase"_orig_file.tsv \
    | tail -n +2 \
    >> "$infilebase"_fulldata.tsv
  rm "$infilebase"_orig_file.tsv
done

for countfile in "$baseout"/*/results/"$sampid"*"$gtype"_fulldata.tsv; do
  rsync "$countfile" "$infilebase"_orig_file.tsv
  awk '{print $5"\t"$3}' "$infilebase"_orig_file.tsv \
    | tail -n +2 \
    >> "$infilebase"_fulldata_barcodes.tsv
  rm "$infilebase"_orig_file.tsv
done

python3 "$script"/combine_samps/sorted_stat_combine_rep.py \
  "$infilebase" \
  "$outfilebase" \
  "$sampid"
  
rsync "$outfilebase"* "$output"/
rm "$infilebase"*
