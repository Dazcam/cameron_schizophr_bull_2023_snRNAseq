import scdrs
import pandas as pd
import scanpy as sc
from anndata import AnnData
import os, subprocess, tempfile
import seaborn as sn

adata = sc.read(snakemake.input[0])
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)


df = pd.DataFrame(adata.obs['n_genes_by_counts'])
df.columns = ['n_genes']
df.index.name = 'index'
df.reset_index(inplace=True)
df.to_csv(snakemake.output[0], sep="\t", index = False)
