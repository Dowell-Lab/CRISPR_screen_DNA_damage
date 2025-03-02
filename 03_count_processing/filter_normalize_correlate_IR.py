#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
from sys import argv
from scipy import stats

### Define inputs/outputs
countfile_ctrl = '/home/lsanford/crispr_screen/counts/fullbasal_grna1_coupled_counts.tsv'
countfile = '/home/lsanford/crispr_screen/counts/IR_grna1_coupled_counts.tsv'
outbase = '/home/lsanford/crispr_screen/finalout/IR'

# countfile_ctrl = '/path/to/coupled/countfile/fullbasal_grna1_coupled_counts.tsv'
# countfile = '/path/to/coupled/countfile/IR_grna1_coupled_counts.tsv'
# outbase = '/path/to/output/directory/<output_basename>'

### Define global variables
# Sample names in countfile header
# Day14 samples should immediately follow Day0 samples of the same replicate
sample_names_to_analyze = [
  'day0_rep1', 'IR_day15_rep1', 'day0_rep2', 'IR_day15_rep2',
] # Omit replicate 3
control_samples = ['day0_rep1', 'day0_rep2',]
replicates = ['rep1', 'rep2',]

### Define functions

def read_file(filename):
  """Read in countfile.

    Parameters:
      filename (str): filename of count file
                      *first line should be header with sample names

    Returns:
      int_df (df):    dataframe of integer counts
  """
  readin_df = pd.DataFrame()
  readin_df = pd.read_csv(
    filename,
    sep = '\t',
    header = 0,
    index_col = 0,
  )
  int_df = readin_df.astype(int) # Convert to integers

  return int_df


def split_guides(df):
  """Split guide identifiers into genes and guides for filtering/grouping.

    Parameters:
      df (df):  dataframe of counts/fold changes/etc.
                *index must be in format <guide1>;<guide2>

    Returns:
      df (df):  data frame with the same data but new columns containing gene
                and guide strs
  """
  df['guides'] = df.index
  df['guide1'] = df.guides.apply(lambda x : x.split(';')[0])
  df['guide2'] = df.guides.apply(lambda x : x.split(';')[1])
  df['gene1'] = df.guide1.apply(lambda x : x.split('_')[0])
  df['gene2'] = df.guide2.apply(lambda x : x.split('_')[0])
  df = df.drop(columns=['guides'])

  return df


def coverage_filter(count_df, mapnum):
  """Filter and normalize counts.

    Parameters:
      count_df (df):        dataframe with raw count data for each replicate
      mapnum (list):        list of proportions for normalization based on 
                            total observed pairs

    Returns:
      norm_count_df (df):   dataframe of filtered and normalized counts
  """
  # Filter out pairs with < 10
  count_keep = count_df[control_samples].apply(
    lambda row: any([(x < 10) for x in row]), axis = 1
  )
  count_keep = ~count_keep
  count_keep = count_keep[count_keep]
  count_df = count_df.loc[count_keep.index]

  # Add pseudocount of 1
  count_df = count_df + 1

  # Normalize to total pairs observed in each sample
  norm_count_df = count_df.div(mapnum, axis='columns')

  return norm_count_df


def fc_calc(norm_count_df):
  """Calculate fold changes.

    Parameters:
      norm_count_df (df): dataframe with normalized count data for
                          each replicate
                          * Assumes day14 follows corresponding day0 sample

    Returns:
      enrich_df (df):     dataframe of fold changes
  """
  fc_df = pd.DataFrame()
  colbound = len(replicates)
  for i in range(0,colbound):
    fc_df[replicates[i]] = norm_count_df.iloc[:, (i * 2 + 1)]/norm_count_df.iloc[:, (i * 2)]

  return fc_df


def filter_fc_norm(guide_count_df):
  """Filter, normalize, and calculate fold changes.

    Parameters:
      guide_count_df (df):      dataframe with raw count data for each
                                individual sample

    Returns:
      guide_fc_df (df):     dataframe of fold changes
  """
  # Calculate normalization factors
  mapped = guide_count_df.sum(axis=0)
  samp_mapped = mapped/np.max(mapped)

  # Run individual steps
  guide_filt_count_df = coverage_filter(guide_count_df, samp_mapped)
  guide_fc_df = fc_calc(guide_filt_count_df)

  # Output counts
  guide_filt_count_df.to_csv(
    (outbase + '_filt_norm_counts.tsv'),
    sep='\t', header=True, index=True
  )

  # Output fold changes
  guide_fc_df.to_csv(
    (outbase + '_filt_norm_fc.tsv'),
    sep='\t', header=True, index=True
  )

  return guide_fc_df


def flip_columns(df, cola, colb):
  """Flip values for dataframe columns.

    Parameters:
      df (df):      dataframe of values
      cola (str):   first df column to switch
      colb (str):   second df column to switch

    Returns:
      df (df):      dataframe with columns flipped
  """
  df['store'] = df[cola]
  df[cola] = df[colb]
  df[colb] = df['store']
  df.drop(columns=['store'], inplace=True)

  return df


def single_phenotype_calc(fc_df):
  """Calculate single-gene phenotypes against nt guides.

    Parameters:
      fc_df (df):         dataframe of fold changes labeled by genes/guides

    Returns:
      single_phen (df):   dataframe of single gene phenotypes
  """
  # Extract all gene guide/nt guide pairs
  # Put all non-nt gene ids in 'gene1' column for grouping
  single_gene_df1 = fc_df.loc[
    (fc_df['gene1'] != 'nt') & (fc_df['gene2'] == 'nt')
  ]
  single_gene_df2 = fc_df.loc[
    (fc_df['gene2'] != 'nt') & (fc_df['gene1'] == 'nt')
  ]
  single_gene_df2 = flip_columns(single_gene_df2, 'gene1', 'gene2')
  single_gene_df2 = flip_columns(single_gene_df2, 'guide1', 'guide2')
  single_gene_df = pd.concat([single_gene_df1, single_gene_df2])

  # Grab list of all genes in df and extract fc values by gene
  genes = single_gene_df['gene1'].unique()
  guides = single_gene_df['guide1'].unique()
  values = single_gene_df.groupby(['gene1', 'guide1'])[replicates]

  # Iterate through genes
  single_gene_phen = pd.Series()
  single_guide_phen = pd.Series()
  omit_guides = []
  for gene in genes:
    all_guides = []
    single_guide = {}
    # Store values for all guide-nt pairs within the gene
    for guide_num in ['_sg1', '_sg2', '_sg3', '_sg4']:
      guide = gene + guide_num
      if guide not in guides:
        continue
      vals = values.get_group((gene, guide))
      single_guide[guide] = list(vals.melt()['value'])
      all_guides = all_guides + list(vals.melt()['value'])
    keep_guides = single_guide.copy()
    # If the fcs for guide-nt pairs are dramatically different
    # for one guide compared to the others for that genes, remove them
    # Otherwise average all fcs for that guide paired with nts and store
    for key, val in single_guide.items():
      pval = stats.ttest_ind(val, all_guides)[1]
      if pval < 0.1:
        omit_guides.append(key)
        keep_guides.pop(key)
      else:
        keep_guides[key] = np.mean(val)
    # Average phenotypes for guide and store as gene-nt phenotype
    # and store guide-nt phenotypes
    single_gene_phen.loc[gene] = np.mean(list(keep_guides.values()))
    for key, val in keep_guides.items():
      single_guide_phen.loc[key] = np.mean(val)

  return single_gene_phen, single_guide_phen


def pair_phenotype_calc(fc_df):
  """Calculate gene pair phenotypes.

    Parameters:
      fc_df (df):         dataframe of fold changes labeled by genes/guides

    Returns:
      paired_phen (df):   dataframe of gene-gene phenotypes
  """
  # Filter to only gene-gene pairs
  guidepair_df = fc_df.loc[
    ~((fc_df['gene1'] == 'nt') | (fc_df['gene2'] == 'nt'))
  ]

  # Average replicates and store
  paired_phen = pd.DataFrame()
  paired_phen['obs_phen'] = guidepair_df[replicates].mean(axis=1)

  return paired_phen


def obs_pred_corrs(obs_pred_df, compkey):
  """Calculate Pearson and Spearman correlations between phenotypes.

    Parameters:
      obs_pred_df (df):  dataframe containing observed and predicted
                         fold changes
      compkey (str):     string of key to specify what to correlate
                         against observed phenotype (should be either
                         'pred_phen_guide' or 'pred_phen_gene')

    Returns:
      statcalc (list):   list of stats as follows:
                         [pearson coeff, pearson pval,
                          spearman coeff, spearman pval]
  """
  stats_df = pd.DataFrame()
  stats_df['logobs'] = np.log2(obs_pred_df['obs_phen'])
  stats_df['logcomp'] = np.log2(obs_pred_df[compkey])
  stats_df = stats_df.loc[
    ~(np.isnan(stats_df['logobs']) | np.isinf(stats_df['logobs']))
  ]
  stats_df = stats_df.loc[
    ~(np.isnan(stats_df['logcomp']) | np.isinf(stats_df['logcomp']))
  ]
  statcalc = list(
      stats.pearsonr(
        x=stats_df['logobs'],
        y=stats_df['logcomp']
      )
    ) + list(
      stats.spearmanr(
        stats_df['logobs'],
        stats_df['logcomp']
      )
    )
  
  return statcalc


def pred_phen_calc(obs_pred_df, single_gene, single_guide):
  """Calculate predicted phenotypes.

    Parameters:
      obs_pred_df (df):  dataframe containing observed fold changes

    Returns:
      obs_pred_df (df):  dataframe with added single-gene/guide fold changes
                    and associated predicted phenotypes
  """
  # Add in single gene and guide phenotypes
  obs_pred_df['gene1_phen'] = obs_pred_df.gene1.apply(
    lambda x : None if x not in single_gene.index else single_gene[x]
  )
  obs_pred_df['guide1_phen'] = obs_pred_df.guide1.apply(
    lambda x : None if x not in single_guide.index else single_guide[x]
  )
  obs_pred_df['gene2_phen'] = obs_pred_df.gene2.apply(
    lambda x : None if x not in single_gene.index else single_gene[x]
  )
  obs_pred_df['guide2_phen'] = obs_pred_df.guide2.apply(
    lambda x : None if x not in single_guide.index else single_guide[x]
  )

  obs_pred_df = obs_pred_df[obs_pred_df['gene1_phen'].notna()]
  obs_pred_df = obs_pred_df[obs_pred_df['gene2_phen'].notna()]

  # Calculate predicted phenotype based on single gene phenotypes
  obs_pred_df['pred_phen_gene'] = obs_pred_df['gene1_phen'] * obs_pred_df['gene2_phen']
  obs_pred_df['pred_phen_guide'] = obs_pred_df['guide1_phen'] * obs_pred_df['guide2_phen']
  obs_pred_df = obs_pred_df.loc[~np.isnan(obs_pred_df['pred_phen_guide'])]
  
  return obs_pred_df


def observed_predicted(single_gene_phen, single_guide_phen, paired_phen):
  """Calculate observed and predicted gene-gene phenotypes.

    Parameters:
      single_gene_phen (df):  dataframe of gene/nt phenotypes (fc)
      single_guide_phen (df): dataframe of guide/nt phenotypes (fc)
      paired_phen (df):       dataframe of gene/gene phenotypes (fc)

    Returns:
      accum_phenotypes (df):  dataframe of observed v predicted gene pair 
                              phenotypes
  """
  # Get gene info and sort
  paired_phen = split_guides(paired_phen)
  paired_phen.sort_values(by=['gene1', 'gene2'],inplace=True)
  unique_gene1 = list(paired_phen.gene1.unique())

  # Initialize dataframes
  accum_phenotypes = pd.DataFrame(columns=[
    'obs_phen','pred_phen_gene','pred_phen_guide','zscore','pval','corr_pval'
  ])
  correlations_genes = pd.DataFrame(columns=[
    'pearson','pearson_pvalue','spearman','spearman_pvalue'
  ])
  correlations_guides = correlations_genes.copy()
  np.seterr(divide='warn')

  for g1 in unique_gene1:
    # Extract all gene/gene pairs for a specific gene and
    # rearrange columns for sorting
    genephen_accum = pd.DataFrame(columns=[
      'obs_phen', 'pred_phen_gene', 'pred_phen_guide'
    ])
    working = paired_phen.loc[paired_phen['gene1'] == g1]
    working2 = paired_phen.loc[paired_phen['gene2'] == g1]
    working2 = flip_columns(working2, 'gene1', 'gene2')
    working2 = flip_columns(working2, 'guide1', 'guide2')
    gene_calcs = pd.concat([working, working2])

    # Calculate predicted phenotypes and correlation statistics
    gene_calcs = pred_phen_calc(gene_calcs, single_gene_phen, single_guide_phen)
    correlations_guides.loc[g1] = obs_pred_corrs(gene_calcs, 'pred_phen_guide')

    # Loop through all g2 paired genes
    unique_gene2 = list(gene_calcs.gene2.unique())
    for g2 in unique_gene2:
      gene_label = ';'.join([g1, g2])

      # Extract guide/guide pairs for the two genes and take phenotype means
      pair = gene_calcs.loc[gene_calcs['gene2'] == g2]
      phenotype_cols = ['obs_phen', 'pred_phen_gene','pred_phen_guide']
      phen_means = pair[phenotype_cols].mean(axis=0)
      
      # Log2 transform phenotypes and store
      for phen in phenotype_cols:
        phen_means[phen] = np.log2(phen_means[phen])
      genephen_accum.loc[gene_label] = phen_means
    
    # Clean data
    for phen in phenotype_cols:
      genephen_accum = genephen_accum.loc[
        (~np.isinf(genephen_accum[phen])) & (~np.isnan(genephen_accum[phen]))
      ]

    # Calculated stats
    genephen_accum['zscore'] = stats.zscore(genephen_accum['obs_phen'])
    genephen_accum['pval'] = genephen_accum.zscore.apply(lambda x : stats.norm.sf(abs(x)))
    genephen_accum['corr_pval'] = genephen_accum['pval'] * len(genephen_accum)
    genephen_accum['corr_pval'].loc[genephen_accum['corr_pval'] > 1] = 1
    correlations_genes.loc[g1] = obs_pred_corrs(genephen_accum, 'pred_phen_guide')

    # Store gene data
    accum_phenotypes = pd.concat([accum_phenotypes, genephen_accum])

  return accum_phenotypes, correlations_genes, correlations_guides


def phenotype_calc(guide_fc_df):
  """Calculate single, paired, observed, and predicted phenotypes.

    Parameters:
      guide_fc_df (df):         dataframe of fold changes for guide pairs

    Returns:
      genepair_phenotypes (df): dataframe of gene-gene phenotypes and significance calcs
  """
  guide_fc_df = split_guides(guide_fc_df)

  # Calculate phenotypes
  gene_phenotype, guide_phenotype = single_phenotype_calc(guide_fc_df)
  guidepair_phenotype = pair_phenotype_calc(guide_fc_df)

  # Calculate predicted phenotypes and stats
  genepair_phenotypes, corr_genes, corr_guides = observed_predicted(
    gene_phenotype, guide_phenotype, guidepair_phenotype
  )

  # Output correlations
  corr_genes.to_csv(
    (outbase + '_gene_obs_pred_pearsons.tsv'),
    sep='\t', index=True
  )
  corr_guides.to_csv(
    (outbase + '_guide_obs_pred_pearsons.tsv'),
    sep='\t', index=True
  )

  genepair_phenotypes = split_guides(genepair_phenotypes)
  genepair_phenotypes = genepair_phenotypes.drop(columns=['guide1', 'guide2'])
  new_col_order = [
    'gene1', 'gene2', 'obs_phen', 'pred_phen_gene', 'pred_phen_guide',
    'zscore', 'pval', 'corr_pval',
  ]
  genepair_phenotypes = genepair_phenotypes[new_col_order]
  genepair_phenotypes.to_csv(
    (outbase + '_phenotypes.tsv'),
    sep='\t', header=True, index=False
  )

  return genepair_phenotypes


def significance_filtering(phenotypes):
  """Calculate single, paired, observed, and predicted phenotypes.

    Parameters:
      phenotypes (df):    dataframe of fold changes for guide pairs

    Returns:
      None
  """
  phenotypes = phenotypes.drop(columns=['zscore', 'corr_pval'])
  phenotypes_copy = phenotypes.copy()

  # Pull reciprocal label and pvalue (for flipped gene pair)
  phenotypes['recip_pair'] = phenotypes.apply(
    lambda x : ';'.join([x.gene2, x.gene1]),
    axis=1
  )
  phenotypes['pval_recip'] = phenotypes.apply(
    lambda x : np.nan if x.recip_pair not in phenotypes_copy.index else phenotypes_copy.loc[x.recip_pair]['pval'],
    axis=1
  )
  phenotypes = phenotypes.loc[~np.isnan(phenotypes['pval_recip'])]

  # Multiply pvalues for each gene in the gene pair and deduplicate df
  phenotypes['comb_pval'] = phenotypes['pval'] * phenotypes['pval_recip']
  non_dup = phenotypes.loc[phenotypes['gene1'] != phenotypes['gene2']]
  pairs = list(non_dup.index)
  pairs.reverse()
  pairs_recip = phenotypes['recip_pair'].to_dict()
  for pair in pairs:
    if pair in list(pairs_recip.values()):
      pairs_recip.pop(pair)
  phenotypes = phenotypes.loc[pairs_recip.keys()]
  phenotypes = phenotypes.drop(columns=['recip_pair'])

  # Output data
  phenotypes.to_csv(
    (outbase + '_reciprocal_phenotypes.tsv'),
    sep='\t', header=True, index=False
  )

### Run script

# Basal
# Read in files

## Read in files
basal_df = read_file(countfile_ctrl)
basal_df = basal_df[control_samples]
IR_df = read_file(countfile)
allcount_df = pd.concat([basal_df, IR_df], axis=1)
allcount_df = allcount_df[sample_names_to_analyze]
allcount_df.replace(np.nan, 0, inplace=True) # Replace any NaNs

## Run pipeline
guide_fold_changes_IR = filter_fc_norm(allcount_df)
phenotypes_IR = phenotype_calc(guide_fold_changes_IR)
significance_filtering(phenotypes_IR)
