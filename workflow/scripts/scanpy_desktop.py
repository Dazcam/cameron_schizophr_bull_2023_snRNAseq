import scdrs
import pandas as pd
import scanpy as sc
from anndata import AnnData
import os, subprocess, tempfile
import seaborn as sn

RES_DIR = '/Users/darren/Desktop/fetal_brain_snRNAseq_GE_270922/resources/raw_data/shi_et_al_2021/'
DATA_DIR = '/Users/darren/Desktop/fetal_brain_snRNAseq_GE_270922/results/'
H5AD_DIR = DATA_DIR + 'h5ad_objects/'
OUT_DIR = '../Desktop/fetal_brain_snRNAseq_GE_270922/results/'

# Load data this was generated in: 
adata = sc.read(H5AD_DIR + 'shi2021_filt.h5ad')

# Access cell / gene names
adata.obs_names
adata.var_names

# Add QC metrics
adata.var['mt'] = adata.var_names.str.startswith('MT-') 
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, 
                           inplace=True)


# Create covariate file
df = pd.DataFrame(adata.obs['n_genes_by_counts'])
df.columns = ['n_genes']
df.index.name = 'index'
df.reset_index(inplace=True)
df.to_csv(OUT_DIR + test.tsv, sep="\t", index = False)


# Trouble shooting 
sc.pp.normalize_per_cell(adata2, counts_per_cell_after=1e4)

adata2 = sc.AnnData(adata)

adata2.write('shi2021_filt.h5ad')




adata3 = sc.read_csv(RES_DIR + 'GSE135827_GE_mat_raw_count_with_week_info.txt', 
                     delimiter= '\t')
adata4 = adata3.transpose()
adata5 = adata

# Add QC metrics
adata4.var['mt'] = adata.var_names.str.startswith('MT-') 
sc.pp.calculate_qc_metrics(adata4, qc_vars=['mt'], percent_top=None, log1p=False, 
                           inplace=True)

