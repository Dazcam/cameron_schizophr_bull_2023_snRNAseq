# -------------------------------------------------------------------------------------
#
#    snRNAseq - scDRS visualisation
#
# -------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Code for figures: clusters lvl 1 and lvl 2
#  Issues: Region IDs encoded differently in MAGMA and LSDR analysis
#          Much tweaking needed between ggplot and cowplot for FDR < 5% line

##  Load Packages  --------------------------------------------------------------------
if (!require("Require")) install.packages("Require")
Require::Require(c("tidyverse", "readxl", "data.table", "ggdendro", "Seurat", 
                   "SeuratWrappers", "cowplot", "scCustomize", "rmarkdown", "SeuratDisk")) 
# scCustomize for contrast in cols


##  Set Variables  --------------------------------------------------------------------
RAW_DIR = '/Users/darren/Desktop/fetal_brain_snRNAseq_GE_270922/resources/raw_data/shi_et_al_2021/'
DATA_DIR <- '/Users/darren/Desktop/fetal_brain_snRNAseq_GE_270922/results/'
R_DIR = paste0(DATA_DIR, 'R_objects/')
scDRS_DIR <- paste0(DATA_DIR, 'scDRS/')

seurat.shi.bc <- readRDS(paste0(R_DIR, 'seurat_shi_bc.rds'))
df_gs <- read_tsv(paste0(scDRS_DIR,'scDRS_genewise_Z_top1K.gs'))
score <- read_tsv(paste0(scDRS_DIR, 'SCZ.full_score.gz'))
level1 <- read_tsv(paste0(scDRS_DIR, 'SCZ.scdrs_group.cluster_level_1'))

score[1:10, 1:10]

seurat.shi.bc$scDRS_score <- score$norm_score
FeaturePlot_scCustom(seurat.shi.bc, na_cutoff = -5, feature = 'scDRS_score', 
                     colors_use = rev(brewer.pal(9, "RdBu"))) 
library(RColorBrewer)

display.brewer.pal(n = 8, name = 'Dark2')
brewer.pal(n = 11, name = "RdBu")
