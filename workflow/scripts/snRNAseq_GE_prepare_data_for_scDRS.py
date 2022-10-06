#------------------------------------------------------------------------------
#
#    Prep Shi 2021 snRNAseq data for scDRS - create h5
#
#------------------------------------------------------------------------------

import scdrs
import pandas as pd
import scanpy as sc
from anndata import AnnData
import os, subprocess, tempfile
import seaborn as sn
from collections import Counter

SHI_GeX = snakemake.input[0]
SHI_META = snakemake.input[1]
H5AD_OUT = snakemake.output[1]
COV_OUT = snakemake.output[0]


# Load data -------------------------------------------------------------------
exp_mat = pd.read_csv(SHI_GeX, delimiter = '\t')

# Metadata
meta = pd.read_excel(SHI_META, header = 1)


##  Check if the cell orders are identical in shi_meta and shi_data  ---------
exp_cellIDs = pd.DataFrame(index = exp_mat.columns)
exp_cellIDs['Cells'] = exp_cellIDs.index
exp_cellIDs[['CellID', 'pcw']] = exp_cellIDs['Cells'].str.split('-GW', expand = True)
exp_cellIDs = exp_cellIDs.reset_index()

# Get the metadata cell IDs and remove trailing number 
meta[['CellID', 'Number']] = meta['Cells'].str.split('-', expand=True)

# Note that the cell IDs in the metadata table do not match that give in the 
# colnames of the data matrix. The former have 1-11 assigned non uniformly
meta.Number.unique()
Counter(meta.Number) # Note these numbers are not uniform

# This number might be important - check if there are any duplicated cell IDs 
# in meta when that number is removed
meta.CellID.duplicated().sum()

# It's likely the numbers refer to different sequencing runs. Cells from 
# different sequencing runs can be given the same ID as the barcodes
# in each GEM are reused. We need to make sure the cell ID for each cell is unique
# Now check cell IDs in shi_data
exp_cellIDs.CellID.duplicated().sum()

# Good sign that duplicate cell counts match - now check if the cell IDs 
# in shi_meta and shi_data are in the same order
exp_cellIDs.CellID.equals(meta.CellID)

# Get metadata into correct format
meta = meta.join(exp_cellIDs.pcw)
meta = meta.set_index("Cells")
meta = meta[["Major types", "pcw"]]
meta.columns = ["ClusterID", "pcw"]

# Some weird behaviour with spaces before progenitor
meta = meta.replace({' progenitor': 'progenitor', '  progenitor': 'progenitor'})

# We've already extracted the pcw info and we know the cell orderings are 
# identical so we can change the cell names in the original shi_data geX 
# matrix to exactly match that in the meta data
exp_mat_t = exp_mat.transpose()
exp_mat_t.index = meta.index


# Convert to scanpy anndata object
adata = AnnData(exp_mat_t, obs = meta)

# Reatin only clusters of interest - 47101 retained (same as in R script)
adata[adata.obs['ClusterID'].isin(['MGE', 'LGE', 'CGE', 'OPC', 'progenitor', 
                                   'Microglia', 'Endothelial'])]

# Access cell / gene names
adata.obs_names
adata.var_names

# Add QC metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-') 
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, 
                           inplace=True)

##  Create covariate file  ----------------------------------------------------
cov_df = pd.DataFrame(adata.obs['n_genes_by_counts'])
cov_df.columns = ['n_genes']
cov_df.index.name = 'index'
cov_df.reset_index(inplace = True)
cov_df.to_csv(COV_OUT, sep = "\t", index = False)

# Test for failure of compute score on hawk
sc.pp.normalize_per_cell(adata, counts_per_cell_after = 1e4)

# Write h5ad
adata.write(H5AD_OUT)

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------