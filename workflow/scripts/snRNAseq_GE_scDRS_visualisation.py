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
score = pd.read_csv(scDRS_DIR + "SCZ.full_score.gz", sep = "\t", index_col = 0)
adata = sc.read(H5AD_DIR + 'shi.bc.qc.h5ad')

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

sc.pl.umap(
    adata, 
    color=["cluster_level_1", "SCZ"],   
    ncols=1,
    color_map="RdBu_r",
    vmin=-5,
    vmax=5
)

sc.pl.umap(adata,
    color='SCZ',
    color_map="RdBu_r",
    vmin=-5,
    vmax=5,
    s=20,
)


## Group level stats ----------------------------------------------------------
dict_df_stats = {
    trait: pd.read_csv(scDRS_DIR + "SCZ.scdrs_group.cluster_level_1", sep="\t", index_col=0)
    for trait in ["SCZ"]
}

dict_celltype_display_name = {
    "MGE": "MGE",
    "LGE": "LGE",
    "CGE": "CGE",
    "Early_InN": "Early_InN",
    "Microglia": "Microglia",
    "Progenitor": "Progenitor",
}

scdrs.util.plot_group_stats(
    {
        trait: df_stats.rename(index = dict_celltype_display_name)
        for trait, df_stats in dict_df_stats.items()
    }
)