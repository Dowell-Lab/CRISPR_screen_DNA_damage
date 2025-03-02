#!/bin/bash

samplist=/path/to/processing/files/samplist_combine.txt
baseout=/path/to/read/count/output/parent/dir
script=/path/to/scripts

# Run sbatch script for each sample
while read line; do
  sbatch --export=sampid=$line,baseout=$baseout "$script"/combine_samps/sorted_stat_combine_rep.sbatch
done < "$samplist"
