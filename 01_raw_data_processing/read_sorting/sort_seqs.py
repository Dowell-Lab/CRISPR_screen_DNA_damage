#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
from sys import argv
import math

## Define inputs/outputs
indir = argv[1]
filebase = argv[2]
outdir = argv[3]
nt_barcode_file = argv[4]

## Define functions
def read_infofile(filename, keyname, exceptions):
  readin = {}
  with open(filename) as f:
    for line in f:
      inline = line.strip().split()
      gene = inline[1]
      if "ZS" in inline[3]:
        first_align = inline[2].split(':')
        second_align = inline[3].split(':')
        if first_align[2] == second_align[2]:
          if gene in exceptions.keys() or gene in exceptions.values():
            if inline[0] in readin.keys():
              continue
            elif gene in exceptions.values():
              readin[inline[0]] = {keyname:gene}
            else:
              readin[inline[0]] = {keyname:exceptions[gene]}
          else:
            continue
        else:
          readin[inline[0]] = {keyname:gene}
      else:
        readin[inline[0]] = {keyname:gene}
  return readin

def read_nt(filename):
  readin = {"seq1":{},"seq2":{}}
  with open(filename) as f:
    for line in f:
      inline = line.strip().split()
      ntnum = int(inline[0].split('_')[1])
      if ntnum < 26:
        readin["seq1"][inline[1]] = inline[0]
      else:
        readin["seq2"][inline[1]] = inline[0]
  return readin

def add_df(readin,full_data):
  adddf = pd.DataFrame(list(readin.values()), index = list(readin.keys()))
  del readin
  full_data = pd.concat([full_data, adddf],axis=1)
  del adddf
  return full_data

def find_match(geneid,seqnum,matching,ntbc):
  if type(geneid) == float and np.isnan(geneid):
    return None
  gspl = geneid.split('_')
  if gspl[0] != "nt":
    if gspl[1] in matching[seqnum].keys():
      gene = gspl[0] + "_" + matching[seqnum][gspl[1]]
    else:
      return None
  elif geneid in ntbc[seqnum].keys():
    gene = ntbc[seqnum][geneid]
  else:
    return None
  return gene

## Define static dicts
guide_except = {
  "Tex15_sg1":"Tex12_sg1",
  "Tex15_sg2":"Tex12_sg2",
  "Tex15_sg3":"Tex12_sg3",
  "Tex15_sg4":"Tex12_sg4",
  "Prkdc_sg1":"Prkcg_sg1",
  "Prkdc_sg2":"Prkcg_sg2",
  "Prkdc_sg3":"Prkcg_sg3",
  "Prkdc_sg4":"Prkcg_sg4",
  "Rmi2_sg2":"Rmi1_sg1",
  "Ino80b_sg4":"Hmga2_sg4",
  "Pold2_sg3":"Pold2_sg2",
}

match_replace = {
  "seq1":{"sg1":"sg2", "sg2":"sg1"},
  "seq2":{"sg3":"sg4", "sg4":"sg3"},
}

nt_barcodes = read_nt(nt_barcode_file)

## Read in files
read_data = read_infofile(
  (indir + filebase + "_grna1.readinfo"),
  "r1_grna",
  guide_except
)
parsed_data = pd.DataFrame(
  list(read_data.values()),
  index = list(read_data.keys())
)
del read_data

read_data = read_infofile(
  (indir + filebase + "_grna2.readinfo"),
  "r2_grna",
  guide_except
)
parsed_data = add_df(read_data,parsed_data)

read_data = read_infofile(
  (indir + filebase + "_barcode1.readinfo"),
  "barcode1",
  guide_except
)
parsed_data = add_df(read_data,parsed_data)

read_data = read_infofile(
  (indir + filebase + "_barcode2.readinfo"),
  "barcode2",
  guide_except
)
parsed_data = add_df(read_data,parsed_data)

# Mark exceptions
parsed_data.loc[
  (parsed_data["r1_grna"] == "Rmi1_sg1") & (parsed_data["barcode1"] == "Rmi2_sg1"),
  "r1_grna"
] = "Rmi2_sg2"
parsed_data.loc[
  (parsed_data["r1_grna"] == "Rmi2_sg1") & (parsed_data["barcode1"] == "Rmi1_sg1"),
  "barcode1"
] = "Rmi2_sg2"
parsed_data.loc[
  (parsed_data["r2_grna"] == "Hmga2_sg4") & (parsed_data["barcode2"] == "Ino80b_sg3"),
  "r2_grna"
] = "Ino80b_sg4"
parsed_data.loc[
  (parsed_data["r2_grna"] == "Ino80b_sg3") & (parsed_data["barcode2"] == "Hmga2_sg4"),
  "barcode2"
] = "Ino80b_sg4"
parsed_data.loc[(parsed_data["r1_grna"] == "Pold2_sg3"), "r1_grna"] = "Pold2_sg2"
parsed_data.loc[(parsed_data["barcode1"] == "Pold2_sg3"), "barcode1"] = "Pold2_sg2"
parsed_data.loc[(parsed_data["r2_grna"] == "Pold2_sg2"), "r2_grna"] = "Pold2_sg3"
parsed_data.loc[(parsed_data["barcode2"] == "Pold2_sg2"), "barcode2"] = "Pold2_sg3"

# Find matches
parsed_data["match1"] = parsed_data.barcode1.apply(
  lambda x : find_match(x, "seq1", match_replace, nt_barcodes)
)
parsed_data["match2"] = parsed_data.barcode2.apply(
  lambda x : find_match(x, "seq2", match_replace, nt_barcodes)
)

# Sort out guide1 only
grna1_only = parsed_data[["r1_grna", "match2", "barcode1", "match1"]].copy()

# Remove missing data
parsed_data = parsed_data[~parsed_data.isnull().any(axis=1)]
grna1_only = grna1_only[~grna1_only.isnull().any(axis=1)]

# Find decoupled data
parsed_data["decoup"] = 0
parsed_data.loc[parsed_data["r1_grna"] != parsed_data["match1"], "decoup"] = 1
parsed_data.loc[parsed_data["r2_grna"] != parsed_data["match2"], "decoup"] = 1
coupled = parsed_data.loc[parsed_data["decoup"] == 0]
coupled = coupled.drop(columns=["match1", "match2", "decoup"])
parsed_data = parsed_data.drop(columns=["match1", "match2", "decoup"])

grna1_only["decoup"] = 0
grna1_only.loc[grna1_only["r1_grna"] != grna1_only["match1"], "decoup"] = 1
grna1_coupled = grna1_only.loc[grna1_only["decoup"] == 0]
grna1_coupled = grna1_coupled.drop(columns=["match1", "decoup"])
grna1_coupled_copy = grna1_coupled.copy()
grna1_coupled.columns = ["r1_grna", "r2_grna", "barcode1"]
grna1_only = grna1_only.drop(columns=["match1", "decoup"])
grna1_only.columns = ["r1_grna", "r2_grna", "barcode1"]

# Figuring out which of grna1 coupled also have grna2 coupled
grna1_coupled_wgrna2 = grna1_coupled_copy.merge(
  parsed_data["r2_grna"],
  how='left',
  left_index=True,
  right_index=True
)
grna1_coupled_wgrna2 = grna1_coupled_wgrna2[~grna1_coupled_wgrna2.isnull().any(axis=1)]
grna1_coupled_wgrna2.to_csv(
  (outdir + filebase + '_grna1_coupled_grna2.tsv'),
  sep='\t',
  header=True,
  index=True
)

grna1_coupled_wgrna2["decoup"] = 0
grna1_coupled_wgrna2.loc[
  grna1_coupled_wgrna2["r2_grna"] != grna1_coupled_wgrna2["match2"],
  "decoup"
] = 1
grna1_coupled_grna2_coupled = grna1_coupled_wgrna2.loc[grna1_coupled_wgrna2["decoup"] == 0]
grna1_coupled_grna2_coupled = grna1_coupled_grna2_coupled.drop(columns=["decoup"])
grna1_coupled_grna2_coupled.to_csv(
  (outdir + filebase + '_grna1_coupled_grna2_coupled.tsv'),
  sep='\t',
  header=True,
  index=True
)

# Write out full sets
parsed_data.to_csv(
  (outdir + filebase + '_g1g2_fulldata.tsv'),
  sep='\t',
  header=True,
  index=True
)

count_table = parsed_data.groupby(["r1_grna", "r2_grna"]).size().reset_index()
count_table.columns = ["guide1", "guide2", "count"]
count_table = count_table.sort_values(by=["count"], ascending=False)
count_table.to_csv(
  (outdir + filebase + '_g1g2_fulldata_guide_count.tsv'),
  sep='\t',
  header=True,
  index=False
)
del parsed_data

grna1_only.to_csv(
  (outdir + filebase + '_grna1_fulldata.tsv'),
  sep='\t',
  header=True,
  index=True
)

count_table = grna1_only.groupby(["r1_grna", "r2_grna"]).size().reset_index()
count_table.columns = ["guide1", "guide2", "count"]
count_table = count_table.sort_values(by=["count"], ascending=False)
count_table.to_csv(
  (outdir + filebase + '_grna1_fulldata_guide_count.tsv'),
  sep='\t',
  header=True,
  index=False
)

#count_table = grna1_only.groupby(["r1_grna", "match1"]).size().reset_index()
#count_table.columns = ["guide1", "match1", "count"]
#count_table = count_table.sort_values(by=["count"], ascending=False)
#count_table.to_csv(
#  (outdir + filebase + '_grna1_fulldata_guidebarcode_count.tsv'),
#  sep='\t',
#  header=True,
#  index=False
#)

del grna1_only

# Write out coupled sets
coupled.to_csv(
  (outdir + filebase + '_g1g2_coupled.tsv'),
  sep='\t',
  header=True,
  index=True
)

count_table = coupled.groupby(["r1_grna", "r2_grna"]).size().reset_index()
count_table.columns = ["guide1", "guide2", "count"]
count_table = count_table.sort_values(by=["count"], ascending=False)
count_table.to_csv(
  (outdir + filebase + '_g1g2_coupled_guide_count.tsv'),
  sep='\t',
  header=True,
  index=False
)
del coupled

grna1_coupled.to_csv(
  (outdir + filebase + '_grna1_coupled.tsv'),
  sep='\t',
  header=True,
  index=True
)

count_table = grna1_coupled.groupby(["r1_grna", "r2_grna"]).size().reset_index()
count_table.columns = ["guide1", "guide2", "count"]
count_table = count_table.sort_values(by=["count"], ascending=False)
count_table.to_csv(
  (outdir + filebase + '_grna1_coupled_guide_count.tsv'),
  sep='\t',
  header=True,
  index=False
)
