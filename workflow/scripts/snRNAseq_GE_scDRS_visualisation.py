import scdrs
import scanpy as sc
sc.set_figure_params(dpi=125)
from anndata import AnnData
from scipy import stats
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
import warnings

RAW_DIR = '/Users/darren/Desktop/fetal_brain_snRNAseq_GE_270922/resources/raw_data/shi_et_al_2021/'
DATA_DIR = '/Users/darren/Desktop/fetal_brain_snRNAseq_GE_270922/results/'
SHI_DIR = RAW_DIR + 'science.abj6641_tables_s2_to_s9/'
scDRS_DIR = DATA_DIR + 'scDRS/'
H5AD_DIR = DATA_DIR + 'h5ad_objects/'

df_gs = pd.read_csv(scDRS_DIR + 'scDRS_genewise_Z_top1K.gs', sep = "\t", index_col = 0)
score = pd.read_csv(scDRS_DIR+ "SCZ.full_score.gz", sep = "\t", index_col = 0)
adata = sc.read(H5AD_DIR + 'shi2021_filt.h5ad')
shi_meta = pd.read_excel(SHI_DIR + 'science.abj6641_table_s2.xlsx', index_col = 0, header = 1)


# Add UMAP to adata
adata.obsm['X_umap'] = shi_meta['UMAP-X']
adata.obsm['Y_umap'] = shi_meta['UMAP-Y']
adata.obsm['Z_umap'] = shi_meta['UMAP-Z']

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)

## Analysis of disease enrichment for individual cells  -----------------------
df_gs = df_gs.loc[
    [
        "SCZ"
    ],
    :,
].rename(
    {
        "SCZ": "SCZ"

    }
)
display(df_gs)

dict_score = {
    'SCZ': pd.read_csv(scDRS_DIR + "SCZ.full_score.gz", sep="\t", index_col=0)

}

for trait in dict_score:
    adata.obs[trait] = dict_score[trait]["norm_score"]

adata.obs['SCZ'] = score["norm_score"]

sc.set_figure_params(figsize=[2.5, 2.5], dpi=150)
sc.pl.umap(
    adata, 
    color=["ClusterID", 'SCZ'],
    ncols=1,
    color_map="RdBu_r",
    vmin=-5,
    vmax=5,
)

sc.pl.umap(adata,
    color='SCZ',
    color_map="RdBu_r",
    vmin=-5,
    vmax=5,
    s=20,
)

        
sc.tl.leiden(adata)
adata.obsm

## Group level stats ----------------------------------------------------------

dict_df_stats = {
    trait: pd.read_csv(scDRS_DIR + "SCZ.scdrs_group.ClusterID", sep="\t", index_col=0)
    for trait in ["SCZ"]
}

dict_celltype_display_name = {
    "MGE": "MGE",
    "LGE": "LGE",
    "CGE": "CGE",
    "OPC": "OPC",
    "Endothelial": "Endothelial",
    "Microglia": "Microglia",
    "progenitor": "Progenitor",
}

scdrs.util.plot_group_stats(
    {
        trait: df_stats.rename(index = dict_celltype_display_name)
        for trait, df_stats in dict_df_stats.items()
    }
)