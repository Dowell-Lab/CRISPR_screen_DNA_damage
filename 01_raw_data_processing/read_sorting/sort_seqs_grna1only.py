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

def add_df(readin,full_data):
  adddf = pd.DataFrame(list(readin.values()), index = list(readin.keys()))
  del readin
  full_data = pd.concat([full_data, adddf],axis=1)
  del adddf
  return full_data

def find_match(geneid,matching):
  gspl = geneid.split('_')
  if gspl[0] != "nt":
    gene = gspl[0] + "_" + matching[gspl[1]]
  elif (int(gspl[1]) % 2) == 0:
    gene = gspl[0] + "_" + str(int(gspl[1]) - 1)
  else:
    gene = gspl[0] + "_" + str(int(gspl[1]) + 1)
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
  "sg1":"sg2",
  "sg2":"sg1",
  "sg3":"sg4",
  "sg4":"sg3",
}

## Read in files
read_data = read_infofile((indir + filebase + "_grna1.readinfo"),"r1_grna",guide_except)
parsed_data = pd.DataFrame(list(read_data.values()), index = list(read_data.keys()))
del read_data

read_data = read_infofile((indir + filebase + "_barcode1.readinfo"),"barcode1",guide_except)
parsed_data = add_df(read_data,parsed_data)

read_data = read_infofile((indir + filebase + "_barcode2.readinfo"),"barcode2",guide_except)
parsed_data = add_df(read_data,parsed_data)

# Remove missing data
parsed_data = parsed_data[~parsed_data.isnull().any(axis=1)]

# Mark exceptions
parsed_data.loc[(parsed_data["r1_grna"] == "Rmi1_sg1") & (parsed_data["barcode1"] == "Rmi2_sg1"),"r1_grna"] = "Rmi2_sg2"
parsed_data.loc[(parsed_data["r1_grna"] == "Rmi2_sg1") & (parsed_data["barcode1"] == "Rmi1_sg1"),"barcode1"] = "Rmi2_sg2"
parsed_data.loc[(parsed_data["r1_grna"] == "Pold2_sg3"),"r1_grna"] = "Pold2_sg2"
parsed_data.loc[(parsed_data["barcode1"] == "Pold2_sg3"),"barcode1"] = "Pold2_sg2"
parsed_data.loc[(parsed_data["barcode2"] == "Pold2_sg2"),"barcode2"] = "Pold2_sg3"

# Mark decoupled
parsed_data["match1"] = parsed_data.barcode1.apply(lambda x : find_match(x,match_replace))
parsed_data["match2"] = parsed_data.barcode2.apply(lambda x : find_match(x,match_replace))

parsed_data["decoup"] = 0
parsed_data.loc[parsed_data["r1_grna"] != parsed_data["match1"],"decoup"] = 1
coupled = parsed_data.loc[parsed_data["decoup"] == 0]
del parsed_data
coupled = coupled.drop(columns=["match1","decoup"])

# Write out coupled sets
coupled.to_csv((outdir + filebase + '_grna1_coupled.tsv'), sep='\t', header=True, index=True)

count_table = coupled.groupby(["r1_grna","match2"]).size().reset_index()
count_table.columns = ["guide1","guide2","count"]
count_table = count_table.sort_values(by=["count"], ascending=False)
count_table.to_csv((outdir + filebase + '_grna1_coupled_guide_count.tsv'), sep='\t', header=True, index=False)

coupled["gene1"] = coupled.r1_grna.apply(lambda x : x.split('_')[0])
coupled["gene2"] = coupled.match2.apply(lambda x : x.split('_')[0])
count_table = coupled.groupby(["gene1","gene2"]).size().reset_index()
count_table.columns = ["gene1","gene2","count"]
count_table = count_table.sort_values(by=["count"], ascending=False)
count_table.to_csv((outdir + filebase + '_grna1_coupled_gene_count.tsv'), sep='\t', header=True, index=False)
