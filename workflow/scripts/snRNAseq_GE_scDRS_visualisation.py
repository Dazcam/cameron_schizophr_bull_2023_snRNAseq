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

DATA_DIR = "/Users/darren/Desktop/fetal_brain_snRNAseq_GE_270922/results/"
scDRS_DIR = DATA_DIR + 'scDRS/'
scDRS_SUBDIR = scDRS_DIR + 'scDRS_shi_bc'
H5AD_DIR = DATA_DIR + 'h5ad_objects/'
REGIONS = ["MGE", "LGE", "CGE", "Progenitor"]
GENE_WINDOWS = ["0UP_0DOWN", "10UP_10DOWN", "10UP_35DOWN", "100UP_100DOWN"]
GENE_WINDOWS = ["10UP_10DOWN"]

## Load data  -----------------------------------------------------------------
# Level 1 - main clusters
for WINDOW in GENE_WINDOWS:
    all_score = pd.read_csv(scDRS_SUBDIR + "/" + WINDOW + 
                            "/SCZ.full_score.gz", sep = "\t", index_col = 0)
    all_group = pd.read_csv(scDRS_SUBDIR + "/" + WINDOW + 
                            "/SCZ.scdrs_group.cluster_level_1", 
                            sep="\t", index_col = 0)
    locals()["all_score_" + WINDOW] = all_score
    locals()["all_group_" + WINDOW] = all_group

all_adata = sc.read(H5AD_DIR + 'shi_bc.h5ad')

# Level 2 - subclusters
for REGION in REGIONS:
    adata = sc.read(H5AD_DIR + "shi_bc_" + REGION + ".h5ad")
    locals()[str(REGION) + "_adata" ] = adata

for WINDOW in GENE_WINDOWS:
    for REGION in REGIONS:
        score = pd.read_csv(scDRS_SUBDIR + '_' + REGION + '/' + WINDOW + '/SCZ.full_score.gz', sep = "\t", index_col = 0)
        group = pd.read_csv(scDRS_SUBDIR + '_' + REGION + '/' + WINDOW + 
                            "/SCZ.scdrs_group.cluster_level_2",
                            sep="\t", index_col = 0)
        locals()[str(REGION) + "_score" + WINDOW] = score
        locals()[str(REGION) + "_group" + WINDOW] = group

del adata, group, score, all_score, all_group

# Top genes
df_gs = pd.read_csv(scDRS_DIR + '/scDRS_genewise_Z_top1K.10UP_10DOWN.gs', sep = "\t", index_col = 0)


## Analysis of disease enrichment for individual cells  -----------------------
dict_score = {
    
    trait: pd.read_csv(f"{scDRS_SUBDIR}_{REGION}/{WINDOW}/{trait}.full_score.gz", sep="\t", index_col=0)
    for trait in df_gs.index

}

for trait in dict_score:
    Progenitor_adata.obs[trait] = dict_score[trait]["norm_score"]

adata.obs['SCZ'] = score["norm_score"]

sc.pl.umap(
    Progenitor_adata, 
    color=["cluster_level_1", "SCZ", "HEIGHT"],   
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
    trait: pd.read_csv(f"{scDRS_SUBDIR}/10UP_10DOWN/{trait}.scdrs_group.cluster_level_1", sep="\t", index_col=0)
    for trait in ["SCZ", "HEIGHT"]
    } 

dict_celltype_display_name = {
    "MGE": "MGE",
    "LGE": "LGE",
    "CGE": "CGE",
    "Early_InN": "Early_InN",
    "Microglia": "Microglia",
    "Progenitor": "Progenitor"
    }


scdrs.util.plot_group_stats(
    {
     trait: df_stats.rename(index = dict_celltype_display_name)
     for trait, df_stats in dict_df_stats.items()
     }
    )

for REGION in REGIONS:
    
    for WINDOW in ['10UP_10DOWN']:

        dict_df_stats = {
            trait: pd.read_csv(f"{scDRS_SUBDIR}/{WINDOW}/{trait}.scdrs_group.cluster_level_2", sep="\t", index_col=0)
            for trait in ["SCZ", "HEIGHT"]
            }
        
        if 'MGE' in REGION:
            
            dict_celltype_display_name = {
                "MGE_0": "MGE_0",
                "MGE_1":  "MGE_1",
                "MGE_2":  "MGE_2",
                "MGE_3":  "MGE_3",
                "MGE_4":  "MGE_4",
                "MGE_5":  "MGE_5",
        
            }
            
        elif 'LGE' in REGION:
                
            dict_celltype_display_name = {
                "LGE_0": "LGE_0",
                "LGE_1": "LGE_1",
                "LGE_2": "LGE_2",
                "LGE_3": "LGE_3",
                "LGE_4": "LGE_4",
                "LGE_5": "LGE_5",
                "LGE_6": "LGE_6",
                "LGE_7": "LGE_7",
                }
            
        elif 'CGE' in REGION:
        
            dict_celltype_display_name = {
                "CGE_0": "CGE_0",
                "CGE_1": "CGE_1",
                "GGE_2": "CGE_2",
                "CGE_3": "CGE_3",
        
            }
            
        else:   
             
            dict_celltype_display_name = {
                "Progenitor_0": "Progenitor_0",
                "Progenitor_1": "Progenitor_1",
                "Progenitor_2": "Progenitor_2",
                "Progenitor_3": "Progenitor_3",
                "Progenitor_4": "Progenitor_4",
                "Progenitor_5": "Progenitor_5",
                "Progenitor_6": "Progenitor_6",
                "Progenitor_7": "Progenitor_7",
                "Progenitor_8": "Progenitor_8",
                "Progenitor_9": "Progenitor_9",
                "Progenitor_10": "Progenitor_10",
                }
    
        scdrs.util.plot_group_stats(
            {
                trait: df_stats.rename(index = dict_celltype_display_name)
                for trait, df_stats in dict_df_stats.items()
            }
            )