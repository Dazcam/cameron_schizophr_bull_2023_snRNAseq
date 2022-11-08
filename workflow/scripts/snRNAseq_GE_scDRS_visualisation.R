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
R_DIR <- paste0(DATA_DIR, 'R_objects/')
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
    scDRS_UMAP <- FeaturePlot_scCustom(SEURAT_OBJ, na_cutoff = NA, features = 'scDRS_score', 
                                      colors_use = rev(brewer.pal(9, "RdBu"))) +
      ggtitle(GWAS)
    
    scDRS_GROUP_UMAP <- plot_grid(CLUSTER_UMAP, scDRS_UMAP)
  
    assign(paste0('scDRS_', REGION, '_', GWAS, '_group_umap'), scDRS_GROUP_UMAP, .GlobalEnv)
    
  }

}


for (REGION in c('all', 'MGE', 'LGE', 'CGE', 'Progenitor')) {
  
  SCZ_df <- get(paste0('scDRS_', REGION, '_SCZ_group')) %>%
    mutate(GWAS = rep('SCZ', times = nrow(.)))
  
  HEIGHT_df <- get(paste0('scDRS_', REGION, '_HEIGHT_group')) %>%
    mutate(GWAS = rep('HEIGHT', times = nrow(.)))
    
  # Should we be doing the FDR correction on the DF seprately or together???
  # Need to cross reference - check scDRS code: https://github.com/martinjzhang/scDRS/blob/master/scdrs/util.py
  ALL_DF <- rbind(SCZ_df, HEIGHT_df) %>%
    mutate(fdr_prop = n_fdr_0.1 / n_cell) %>%
    mutate(assoc_mcp_adj = p.adjust(assoc_mcp, method='BH')) %>%
    mutate(assoc_hetero_adj = p.adjust(hetero_mcp, method='BH')) %>%
    mutate(assoc_mcp_bool = ifelse(assoc_mcp_adj < 0.05, TRUE, FALSE)) %>%  
    mutate(hetero_mcp_bool = ifelse(assoc_hetero_adj < 0.05, TRUE, FALSE))
  
  scDRS_GROUP_HEATMAP <- ggplot(data = ALL_DF, mapping = aes(x = group, y = GWAS, fill = fdr_prop, color = 'black')) +
    geom_tile(color = "grey", lwd = 0.7, linetype = 1) +
    geom_point(aes(size = ifelse(assoc_mcp_bool, "square", "no_square")), shape = 0, color = 'black') +
    geom_point(aes(size = ifelse(assoc_mcp_bool, "cross", "no_cross")), shape = 4, color = 'black') +
    scale_size_manual(values = c(square = 10, no_square = NA, cross = 8, no_cross = NA), guide = "none") +
    scale_fill_gradient(name = 'Prop. of sig. cells (FDR < 0.1)', 
                        low = "white", high = "#8B0000", 
                        breaks=c(0, 0.1, 0.2), limits=c(0, 0.2),
                        guide_legend(title.position = "right")) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.ontop = TRUE,
          panel.background = element_rect(fill = "transparent"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          plot.title = element_text(hjust = 0.5, face = 'bold'),
          axis.title.x = element_text(colour = "#000000", size = 14),
          axis.title.y = element_text(colour = "#000000", size = 14),
          axis.text.x  = element_text(colour = "#000000", size = 13, vjust = 0.5, angle = 45),
          axis.text.y  = element_text(colour = "#000000", size = 13),
          axis.ticks = element_blank(),
          legend.position = "top") +
    xlab(NULL) +
    ylab(NULL) +
    coord_fixed()
  
  assign(paste0('scDRS_', REGION, '_group_heatmap'), scDRS_GROUP_HEATMAP, .GlobalEnv)
  assign(paste0('scDRS_', REGION, '_group_corr_df'), ALL_DF, .GlobalEnv)

}

scDRS_group_heatmap <- plot_grid(scDRS_all_group_heatmap, scDRS_MGE_group_heatmap,
          scDRS_CGE_group_heatmap, scDRS_LGE_group_heatmap,
          scDRS_Progenitor_group_heatmap, align = 'tblr')
# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------

