#--------------------------------------------------------------------------------------
#
#    Shi 2021 snRNAseq data - Functions
#
#--------------------------------------------------------------------------------------

## Functions  -------------------------------------------------------------------------

# Test and plots for various resolution parameters
seurat_resolution_test <- function(SEURAT_OBJECT,
                                   CLUSTER_ID,
                                   RESOLUTION,
                                   GROUPBY_VAR,
                                   MARKERS,
                                   UNIQUE_ID) {
  
  for (RES in RESOLUTION) {
    
    # Set the default Ident
    Idents(object = SEURAT_OBJECT) <- SEURAT_OBJECT@meta.data[paste0("RNA_snn_res.", RES)]
    
    # Stats - scCustomize package for stats and nice plots
    pct_stats <- Cluster_Stats_All_Samples(seurat_object = SEURAT_OBJECT, group_by_var = GROUPBY_VAR)
    cluster_stats <- Cluster_Stats_All_Samples(seurat_object = SEURAT_OBJECT, group_by_var = GROUPBY_VAR) %>%
      select_if(!grepl("%", names(.))) 
    median_stats_pcw <- Median_Stats(seurat_object = SEURAT_OBJECT, group_by_var = GROUPBY_VAR)
    median_stats_clust <- Median_Stats(seurat_object = SEURAT_OBJECT, group_by_var = paste0("RNA_snn_res.", RES))
    
    # Plots
    pcw_plot <- DimPlot_scCustom(SEURAT_OBJECT, group.by = GROUPBY_VAR, 
                                 DiscretePalette_scCustomize(num_colors = 26, palette = "ditto_seq"))
    clust_plot <- DimPlot_scCustom(SEURAT_OBJECT, group.by = paste0("RNA_snn_res.", RES),
                                   DiscretePalette_scCustomize(num_colors = 26, palette = "polychrome"))
    shiID_plot <- DimPlot_scCustom(SEURAT_OBJECT, group.by = CLUSTER_ID,
                                   DiscretePalette_scCustomize(num_colors = 26, palette = "polychrome"))
    violin_geX2 <- VlnPlot(SEURAT_OBJECT, MARKERS, stack = TRUE, flip = TRUE, 
                           cols = DiscretePalette_scCustomize(num_colors = 26, palette = "polychrome"),
                           same.y.lims = TRUE, fill.by = 'ident', group.by = paste0("RNA_snn_res.", RES))
    dot_plot <- DotPlot_scCustom(seurat_object = SEURAT_OBJECT, features = rev(MARKERS), 
                                 colors_use = viridis_plasma_dark_high, flip_axes = TRUE, 
                                 group.by = paste0("RNA_snn_res.", RES))
    
    pct_plot <- pct_stats %>%
      select_if(grepl("%|Cluster", names(.))) %>%
      filter(!Cluster %in% 'Total') %>%
      pivot_longer(!Cluster, names_to = GROUPBY_VAR, values_to = "pct") %>%
      mutate(Cluster = factor(Cluster, levels = seq(0, length(unique(Cluster)) - 1, 1))) %>%
      ggplot(aes(fill = pcw, y = pct, x = Cluster)) + 
      geom_bar(position = "fill", stat = "identity") +
      scale_fill_manual(values = DiscretePalette_scCustomize(num_colors = 26, palette = "ditto_seq")) +
      theme_bw()
    
    group_plot <- cowplot::plot_grid(pcw_plot, clust_plot, shiID_plot, 
                                     pct_plot, violin_geX2, dot_plot)
    
    
    # UMAPs
    mge_feat_umap <- FeaturePlot(SEURAT_OBJECT, features = mge_markers)
    lge_feat_umap <- FeaturePlot(SEURAT_OBJECT, features = lge_markers)
    cge_feat_umap <- FeaturePlot(SEURAT_OBJECT, features = cge_markers, ncol = 2)
    left_col <- plot_grid(mge_feat_umap, lge_feat_umap, labels = c('MGE', 'LGE'), 
                          label_size = 12, ncol = 1)
    ge_group_umap <- plot_grid(left_col, cge_feat_umap, labels = c('', 'CGE'), 
                               label_size = 12, ncol = 2)
    
    shi1C_feat_umap <- FeaturePlot(SEURAT_OBJECT, features = shi1C_markers)
    supp1C_feat_umap <- FeaturePlot(SEURAT_OBJECT, features = supp1C_markers)
    our_feat_umap <- FeaturePlot(SEURAT_OBJECT, features = our_markers)
    
    assign(paste0(UNIQUE_ID, '_cluster_stats_', RES), cluster_stats, envir = .GlobalEnv)
    assign(paste0(UNIQUE_ID, '_median_stats_pcw_', RES), median_stats_pcw, envir = .GlobalEnv)
    assign(paste0(UNIQUE_ID, '_median_stats_clust_', RES), median_stats_clust, envir = .GlobalEnv)
    assign(paste0(UNIQUE_ID, '_group_plot_', RES), group_plot, envir = .GlobalEnv)
    
    assign(paste0(UNIQUE_ID, '_ge_group_umap_', RES), ge_group_umap, envir = .GlobalEnv)
    assign(paste0(UNIQUE_ID, '_shi1C_umap_', RES), shi1C_feat_umap, envir = .GlobalEnv)
    assign(paste0(UNIQUE_ID, '_supp1C_umap_', RES), supp1C_feat_umap, envir = .GlobalEnv)
    assign(paste0(UNIQUE_ID, '_our_feat_umap_', RES), our_feat_umap, envir = .GlobalEnv)
    
    if (exists('REGION')) {
    
      if (REGION == 'CGE') {
        
        lvl2_umap <- FeaturePlot(SEURAT_OBJECT, features = cge_markers)
        assign(paste0(UNIQUE_ID, '_lvl2_umap_', RES), lvl2_umap, envir = .GlobalEnv)
        
      } else {
        
        lvl2_umap <- FeaturePlot(SEURAT_OBJECT, features = get(paste0(tolower(REGION), '_level2_markers')))
        assign(paste0(UNIQUE_ID, '_lvl2_umap_', RES), lvl2_umap, envir = .GlobalEnv)
        
      }
    
    } 
    
  }
  
}

## GENE LISTS
# ATAC gene lists
MGE_GENES <- c('CRABP1', 'NKX2-1', 'ETV1', 'NFIA', # Striatal and CRT Ns
               'MEF2C', 'MAF', 'ARX', 'ZEB2', # Cortical interneurons
               'ZFHX3', # Subpallial and GABA Ns also tested 'NR2F1', 'NR2F2', 
               'LHX8', 'CNTNAP2', 'GRIA2', 'ZIC1', # Subpallial Cholin Ns also tested 'ISL1', 'GBX2',
               'MDK', 'LHX6', 'SOX6', 'CXCR4', # Rest from Shi fig 5B
               'ERBB4', 'ANGPT2', 
               'SST'
               
)

LGE_GENES <- c(#'MEIS2', 'TLE4', 'FOXP1', # 'ZNF503' (duplicate), # Striatal potential
  'CHD7', 'DLX5', 'ID2', 'PAX6', # OB precursor
  'PENK', 'CXCL12', 'SP9', # D2 MSN precursor
  'TAC1', 'EBF1', 'ISL1', 'ZNF503', 'PDYN', # D1 MSN precursor (PDYN)
  'TSHZ1', 'FOXP2', 'PBX3', 'ERBB4' #, # D1 MSN precursor (TSHZ1)
  #'SIX3', 'ID4', 'ZFHX4', 'PHLDA1', # Rest from Shi fig 4H
  #'GAP43', 'ZNF521', 'NTRK3',
  #'LHFP', 'ZADH2', 'PBX1',
  #'ADRA2A', 'DRD1', 'DRD2'
)

# LGE_GENES <- c('FOXP1', 'FOXP2', 'ISL1', 'SERTAD4', 'ZNF503',
#                'SIX3', 'ZFHX3', 'SCGN', 'PBX1', 'ZNF521',
#                'PDYN', 'EBF1', 'GRIA2', 'CNTNAP2', 'ZFHX3')

CGE_GENES <- c(
  'NRIP3', 'IER2', 'PRKCA', 'ERBB4',  #CGE-InN-1
  'VSTM2A', 'BEX2', 'PRKCA', 'ERBB4', 'ARL6IP5', # CGE-InN-2
  'ARL4D', 'HIGD2A', 'GNG5', 'SP9', 'STMN4','PBX1', 'MDK', 'SOX6', 'HES6',
  'GAD2', 'CALB2' #CGE-InN-3
)

PROGENITOR_GENES <- c('SOX6', 'DLX6', 'SFTA3', 'NKX2-1', # 'LHX6', 
                      'RIC3', 'MBIP', 'OLIG2', 'PCDH17', # 'GSX1',
                      'NR2F1', 'NR2F2', 'NTRK2', 'DLX5', 'RBP1', 
                      'ZFHX3', 'MEIS2', # 'VIP', 'CALB2', 'SCGN',
                      'SIX3', 'PHLDA1', 'GLI3', 'LGALS1', # 'ZNF503',
                      'PAX6', 'GSX2', 'HEY1'
                      
)

IPC_GENES <- c('DLL1', 'DLL3', 'CCND2', 'GAD1', 'GAD2') # 'KI67', 'EOMES'

LEVEL_1_MARKERS <- c("CALB2", "SCGN", "PCDH9", "ANKS1B",   # CGE
                     "FOXP1", "ZNF503", "SERTAD4", "ISL1", # LGE
                     "LHX6", "NXPH1",                      # MGE
                     "GAD2",                               # InN
                     "HES1", 'OLIG2', "PAX6", "GSX2", 'PTPRZ1', # Progenitor
                     "DLL1", "DLL3", "CCND2",              # IPC
                     "SPI1", "CD68")                       # MG
