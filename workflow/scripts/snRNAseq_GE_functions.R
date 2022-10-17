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
    
  }
  
}