#!/usr/bin/env python
# coding: utf-8

import pandas as pd
from sys import argv
from pathlib import Path

## Define inputs/outputs
indir = argv[1]
outfilebase = argv[2]
suffix = argv[3]

## Define functions
def read_combfile(filename):
  readin = {}
  firstline = 0
  with open(filename) as f:
    for line in f:
      if firstline < 1:
        firstline += 1
        inline = line.strip().split()
        sample = inline[1]
        continue
      else:
        inline = line.strip().split()
        readin[inline[0]] = inline[1]
  return readin,sample

def add_df(readin,samp,full_data):
  adddf = pd.DataFrame({samp:list(readin.values())}, index = list(readin.keys()))
  del readin
  full_data = pd.concat([full_data, adddf],axis=1)
  del adddf
  return full_data

# Find files
if 'IR' in outfilebase:
  pathlist = list(Path(indir).glob('**/IR*' + suffix))
else:
  pathlist = list(Path(indir).glob('**/day0*' + suffix))
  pathlist = pathlist + list(Path(indir).glob('**/day14*' + suffix))

# Read files and combine
filenum = 0
samples = []
for path in pathlist:
  read_data,sampname = read_combfile(str(path))
  samples.append(sampname)
  if filenum < 1:
    filenum += 1
    parsed_data = pd.DataFrame({sampname:list(read_data.values())}, index = list(read_data.keys()))
    del read_data
  else:
    parsed_data = add_df(read_data,sampname,parsed_data)

# Write out full file
samples.sort()
parsed_data = parsed_data[samples]
parsed_data = parsed_data.fillna(0)
parsed_data.to_csv((outfilebase + suffix), sep='\t', header=True, index=True)
