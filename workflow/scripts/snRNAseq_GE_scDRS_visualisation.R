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
                   "SeuratWrappers", "cowplot", "scCustomize", "rmarkdown", "RColorBrewer")) 
# scCustomize for contrast in cols


##  Set Variables  --------------------------------------------------------------------
DATA_DIR <- '/Users/darren/Desktop/fetal_brain_snRNAseq_GE_270922/results/'
R_DIR = paste0(DATA_DIR, 'R_objects/')
scDRS_DIR <- paste0(DATA_DIR, 'scDRS/')
scDRS_SUBDIR <- paste0(scDRS_DIR, 'scDRS_shi_bc_')

##  Load scDRS objects   --------------------------------------------------------------
df_gs <- read_tsv(paste0(scDRS_DIR,'scDRS_genewise_Z_top1K.10UP_10DOWN.gs'))

for (REGION in c('', 'MGE', 'LGE', 'CGE', 'Progenitor')) {
  
  for (GWAS in c('SCZ', 'HEIGHT')) {

  cat('\nLoading scDRS objects for:', REGION, GWAS, '...\n\n')

    if (REGION == '') {
      
      SEURAT_OBJ <- readRDS(paste0(R_DIR, 'seurat_shi_bc', REGION, '.rds'))
      SCORES <- read_tsv(paste0(scDRS_DIR, 'scDRS_shi_bc/10UP_10DOWN/', GWAS, '.full_score.gz'))
      CLUSTS <- read_tsv(paste0(scDRS_DIR, 'scDRS_shi_bc/10UP_10DOWN/', GWAS, '.scdrs_group.cluster_level_1'))
      
      # Add normalised scDRS scores to Seurat object
      SEURAT_OBJ$scDRS_score <- SCORES$norm_score
      
      assign(paste0('scDRS_all_', GWAS, '_seurat'), SEURAT_OBJ, .GlobalEnv)
      assign(paste0('scDRS_all_', GWAS, '_scores'), SCORES, .GlobalEnv)
      assign(paste0('scDRS_all_', GWAS, '_group'), CLUSTS, .GlobalEnv)
      
    } else {
      
      SEURAT_OBJ <- readRDS(paste0(R_DIR, 'seurat_shi_bc_', REGION, '.rds'))
      SCORES <- read_tsv(paste0(scDRS_SUBDIR, REGION, '/10UP_10DOWN/', GWAS, '.full_score.gz'))
      CLUSTS <- read_tsv(paste0(scDRS_SUBDIR, REGION, '/10UP_10DOWN/', GWAS, '.scdrs_group.cluster_level_2'))
      
      # Add scDRS normalised scDRS scores to Seurat objec
      SEURAT_OBJ$scDRS_score <- SCORES$norm_score
      
      assign(paste0('scDRS_', REGION, '_', GWAS, '_seurat'), SEURAT_OBJ, .GlobalEnv)
      assign(paste0('scDRS_', REGION, '_', GWAS, '_scores'), SCORES, .GlobalEnv)
      assign(paste0('scDRS_', REGION, '_', GWAS, '_group'), CLUSTS, .GlobalEnv)
      
    }
    
  }
  
}

for (REGION in c('all', 'MGE', 'LGE', 'CGE', 'Progenitor')) {
  
  for (GWAS in c('SCZ', 'HEIGHT')) {
  
    SEURAT_OBJ <- get(paste0('scDRS_', REGION, '_', GWAS, '_seurat'))
    CLUSTER_UMAP <- DimPlot_scCustom(SEURAT_OBJ, reduction = 'umap', 
                                     group.by = 'cluster_level_2') +
      ggtitle(REGION)
    scDRS_UMAP <- FeaturePlot_scCustom(SEURAT_OBJ, na_cutoff = NA, 
                                      feature = 'scDRS_score', 
                                      colors_use = rev(brewer.pal(9, "RdBu"))) +
      ggtitle(GWAS)
    
    scDRS_GROUP_UMAP <- plot_grid(CLUSTER_UMAP, scDRS_UMAP)
  
    assign(paste0('scDRS_', REGION, '_', GWAS, '_group_umap'), scDRS_GROUP_UMAP, .GlobalEnv)
    
  }

}


# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------

