#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
from sys import argv
import math

## Define inputs/outputs
indir = argv[1]
outdir = argv[2]
sampid = argv[3]

## Define functions
def readwrite_statfile(infile,outfile,colnames):
  indata = pd.read_csv(
    infile,sep='\t',header=None,names=colnames,index_col=False
  )
  collapsed = indata.groupby(colnames).size().reset_index()
  collapsed = collapsed.rename(columns={0:sampid})
  collapsed['guidepair'] = collapsed.apply(
    lambda x: ';'.join([x[colnames[0]],x[colnames[1]]]),axis=1
  )
  collapsed = collapsed.drop(columns = colnames)
  collapsed = collapsed[['guidepair',sampid]]
  collapsed.to_csv(outfile,sep='\t',header=True,index=False)

## Read in files
readwrite_statfile(
  (indir + "_coupled.tsv"),
  (outdir + "_coupled.tsv"),
  ['guide1','bc2_guide']
)
readwrite_statfile(
  (indir + "_fulldata.tsv"),
  (outdir + "_fulldata.tsv"),
  ['guide1','bc2_guide']
)
readwrite_statfile(
  (indir + "_fulldata_barcodes.tsv"),
  (outdir + "_fulldata_barcodes.tsv"),
  ['bc1_guide','bc2_guide']
)
