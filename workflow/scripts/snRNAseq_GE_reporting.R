# -------------------------------------------------------------------------------------
#
#    snRNAseq - Reporting
#
# -------------------------------------------------------------------------------------

##  Set Variables  --------------------------------------------------------------------
SCRIPT_DIR <- '~/Desktop/fetal_brain_snRNAseq_GE_270922/workflow/scripts/'

##  Load objects  --------------------------------------------------------------------
source(paste0(SCRIPT_DIR, 'snRNAseq_GE_plots_magma_and_ldsr_plots.R'))
source(paste0(SCRIPT_DIR, 'snRNAseq_GE_scDRS_visualisation.R'))
source(paste0(SCRIPT_DIR, 'snRNAseq_GE_functions.R'))
source(paste0(SCRIPT_DIR, 'snRNAseq_GE_marker_genes.R'))

# Create markdown file
OUT_DIR <- '~/Desktop/fetal_brain_snRNAseq_GE_270922/results/'
R_DIR <- paste0(OUT_DIR, 'R_objects/')
MARKDOWN_FILE <- paste0(SCRIPT_DIR, 'snRNAseq_GE_reporting.Rmd')
REPORT_DIR <- paste0(OUT_DIR, 'rmarkdown_reports/')
EXCEL_DIR <- paste0(OUT_DIR, 'excel_sheets/')
REPORT_FILE <- 'snRNAseq_GE_reporting.html'
RUN_DIFF_GEX <- FALSE  # Run differential gene expression


##  Load Seurat Objects  --------------------------------------------------------------
seurat.shi.bc <- readRDS(paste0(R_DIR, 'seurat_shi_bc.rds'))
seurat.shi.bc_LGE <- readRDS(paste0(R_DIR, 'seurat_shi_bc_LGE.rds'))
seurat.shi.bc_MGE <- readRDS(paste0(R_DIR, 'seurat_shi_bc_MGE.rds'))
seurat.shi.bc_CGE <- readRDS(paste0(R_DIR, 'seurat_shi_bc_CGE.rds'))
seurat.shi.bc_Progenitor <- readRDS(paste0(R_DIR, 'seurat_shi_bc_Progenitor.rds'))
seurat.shi.bc_dwnSmpl_lvl1 <- readRDS(paste0(R_DIR, 'seurat_shi_bc_dwnSmpl_lvl1.rds'))
seurat.shi.bc_dwnSmpl_lvl2 <- readRDS(paste0(R_DIR, 'seurat_shi_bc_dwnSmpl_lvl2.rds'))


## Subset Shi data - sub-clusters  ---------------------------------------------
for (REGION in c('MGE', 'LGE', 'CGE', 'Progenitor')) {
  
  cat('\nRunning Seurat sub-cluster analysis for:', REGION, '...\n\n')
  
  if (REGION == 'MGE') {
    
    MARKERS <- mge_level2_markers
    
  } else if (REGION == 'LGE') {
    
    MARKERS <- lge_level2_markers
    
  } else if (REGION == 'CGE') {
    
    MARKERS <- cge_markers
    
  } else {
    
    MARKERS <- progenitor_level2_markers 
    
  }
  
  seurat_resolution_test(get(paste0('seurat.shi.bc_', REGION)), 
                         'ClusterID',
                         c(0.5),
                         'pcw',
                         MARKERS,
                         paste0('shi_', REGION))
  
}


## Similarity between Shi level 1 clusters and ours  ----------------------------------
# https://davetang.org/muse/2017/09/28/rand-index-versus-adjusted-rand-index/
shi_IDs <- as.data.frame(seurat.shi.bc$ClusterID)
shi_IDs <- cbind(shi_IDs, seurat.shi.bc$cluster_level_1) %>%
  mutate(across(everything(), ~str_replace(., "MGE", "1" ))) %>%
  mutate(across(everything(), ~str_replace(., "LGE", "2" ))) %>%
  mutate(across(everything(), ~str_replace(., "CGE", "3" ))) %>%
  mutate(across(everything(), ~str_replace(., "Progenitor", "4" ))) %>%
  mutate(across(everything(), ~str_replace(., "Microglia", "5" ))) %>% 
  mutate(across(everything(), ~str_replace(., "Early_InN", "6" ))) %>% 
  mutate(across(everything(), ~str_replace(., "Endothelial", "7" ))) %>%
  mutate(across(everything(), ~str_replace(., "OPC", "8" ))) %>%
  rename(Cluster_ID = `seurat.shi.bc$ClusterID`,
         cluster_level_1 = `seurat.shi.bc$cluster_level_1`)

# Run adjusted rand index for similarity for cluster mappings
truth    <- as.numeric(shi_IDs %>% pull(Cluster_ID))
estimate <- as.numeric(shi_IDs %>% pull(cluster_level_1))
fossil::adj.rand.index(truth, estimate)

##  DeX genes  ------------------------------------------------------------------------
# Get level 2 DeX
if (RUN_DIFF_GEX == TRUE) {

  Idents(seurat.shi.bc) <- seurat.shi.bc$cluster_level_2 
  dEx.markers <- FindAllMarkers(seurat.shi.bc, only.pos = TRUE, 
                                min.pct = 0.25, logfc.threshold = 0.25)
  
  dEx.markers_ordered <- dEx.markers %>%
    arrange(cluster) %>%
    separate(cluster, c('region', 'number')) %>%
    arrange(region, number) %>%
    unite("cluster", region:number) %>%
    mutate_at('cluster', str_replace, "_NA", "") 
  
  dEx.markers_list <- dEx.markers_ordered %>%
    group_split(cluster) %>%
    setNames(unique(dEx.markers_ordered$cluster)) %>%
    writexl::write_xlsx(paste0(EXCEL_DIR, 'snRNAseq_GE_diff_GEx_genes_level_2.xlsx'))
  
  # Replace Idents
  Idents(seurat.shi.bc) <- seurat.shi.bc$cluster_level_1 
  
}

##  Cell type predictions  ------------------------------------------------------------
# MGE - Str-Ctx-InNs - Striatal and cortical interneurons
# MGE - SubP-Chol-InNs - Subpallial cholinergic neurons
# MGE - Ctx-InNs - Cortical interneurons

# LGE - OB-N - Olfactory Bulb neuron precursor
# LGE - D1-MSN-A - D1 MSN (PDYN) precursor
# LGE - D1-MSN-B - D1 MSN (TSHZ1) precursor
# LGE - D2-MSN - D2 MSN precursor

lvl2_cell_IDs <- dput(sort(unique(seurat.shi.bc$cluster_level_2)))
lvl2_cell_labels <- c("CGE_0", "CGE_1", "CGE_2", "CGE_3", "OB-N", "D1-MSN-A", "D2-MSN",
                      "D1-MSN-A", "D2-MSN", "LGE_5", "LGE_6", "D1-MSN-A", "MGE_0", "MGE_1", 
                      "SubP-Chol-InNs", "Ctx-InNs", "MGE_4", "Str-Ctx-InNs", "Other", "Progenitor_0", 
                      "Progenitor_1", "Progenitor_10", "Progenitor_2", "Progenitor_3", 
                      "Progenitor_4", "Progenitor_5", "Progenitor_6", "Progenitor_7", 
                      "Progenitor_8", "Progenitor_9")
lvl2_cell_types <-
  tibble(Original_IDs = lvl2_cell_IDs,
         New_IDs = lvl2_cell_labels)

##  MAGMA & LDSR  ---------------------------------------------------------------------
# Default MAGMA 35UP_10DOWN, LDSR 100UP_100DOWN
magma_ldsr_lvl1_plot <- plot_grid(SCZ_magma_ldsr_lvl_1_plot, HEIGHT_magma_ldsr_lvl_1_plot)
magma_ldsr_lvl2_plot <-plot_grid(SCZ_magma_ldsr_lvl_2_plot, HEIGHT_magma_ldsr_lvl_2_plot)
magma_ldsr_cond_lvl2_plot <- plot_grid(SCZ_magma_ldsr_lvl_2_plot, 
                                       plot_grid(magma_cond_GE_skene_MSN_35UP_10DOWN_plot, 
                                                 magma_cond_GE_skene_InN_35UP_10DOWN_plot, ncol = 1))
magma_ldsr_cond_lvl2_internal_plot <- plot_grid(magma_cond_GE_CGE_1_35UP_10DOWN_plot, magma_cond_GE_CGE_2_35UP_10DOWN_plot,
                                           magma_cond_GE_LGE_0_35UP_10DOWN_plot, magma_cond_GE_LGE_2_35UP_10DOWN_plot,
                                           magma_cond_GE_LGE_4_35UP_10DOWN_plot, magma_cond_GE_MGE_2_35UP_10DOWN_plot,
                                           magma_cond_GE_MGE_3_35UP_10DOWN_plot)

# Correlation between cell number and MAGMA / LDSR results?
# Cluster level 1
cnts_lvl1 <- as.data.frame(table(seurat.shi.bc$cluster_level_1))
cnts_lvl1_dwnSmpl <- as.data.frame(table(seurat.shi.bc_dwnSmpl_lvl1$cluster_level_1))
colnames(cnts_lvl1) <- c('Category', 'cnts')
colnames(cnts_lvl1_dwnSmpl) <- c('Category', 'cnts_ds')
merge_lvl1 <- magma_SCZ_lvl_1_35UP_10DOWN_df %>% 
  left_join(cnts_lvl1) %>%
  left_join(cnts_lvl1_dwnSmpl) %>%
  left_join(ldsr_SCZ_lvl_1_100UP_100DOWN_df) 

cor(merge_lvl1[, unlist(lapply(merge_lvl1, is.numeric))])

cor.test(merge_lvl1$MAGMA, merge_lvl1$cnts)
cor.test(merge_lvl1$LDSR, merge_lvl1$cnts)
plot(merge_lvl1$MAGMA, merge_lvl1$cnts)
plot(merge_lvl1$MAGMA, merge_lvl1$cnts)

# Cluster level 2
cnts_lvl2 <- as.data.frame(table(seurat.shi.bc$cluster_level_2))
cnts_lvl2_dwnSmpl <- as.data.frame(table(seurat.shi.bc_dwnSmpl_lvl2$cluster_level_2))
colnames(cnts_lvl2) <- c('Category', 'cnts')
colnames(cnts_lvl2_dwnSmpl) <- c('Category', 'cnts_ds')
merge_lvl2 <- magma_SCZ_lvl_2_35UP_10DOWN_df %>% 
  left_join(cnts_lvl2) %>%
  left_join(cnts_lvl2_dwnSmpl) %>%
  left_join(ldsr_SCZ_lvl_2_100UP_100DOWN_df) 
  
cor(merge_lvl1[, unlist(lapply(merge_lvl1, is.numeric))])

for (LEVEL in c(1, 2)) {

  MAGMA_PRE_DS <- ggpubr::ggscatter(get(paste0('merge_lvl', LEVEL)), x = "MAGMA", y = "cnts", 
                            add = "reg.line", conf.int = TRUE, 
                            cor.coef = TRUE, cor.method = "pearson",
                            xlab = "MAGMA", ylab = "Cell counts", 
                            title = paste0('Cluster level ', LEVEL))
  
  LDSR_PRE_DS <- ggpubr::ggscatter(get(paste0('merge_lvl', LEVEL)), x = "LDSR", y = "cnts", 
                                   add = "reg.line", conf.int = TRUE, 
                                   cor.coef = TRUE, cor.method = "pearson",
                                   xlab = "LDSR", ylab = "Cell counts", 
                                   title = paste0('Cluster level ', LEVEL))
  
  MAGMA_POST_DS <- ggpubr::ggscatter(get(paste0('merge_lvl', LEVEL)), x = "MAGMA", y = "cnts_ds",
                                     add = "reg.line", conf.int = TRUE, 
                                     cor.coef = TRUE, cor.method = "pearson",
                                     xlab = "MAGMA", ylab = "Cell counts", 
                                     title = paste0('Cluster level ', LEVEL, ' downsampled')) +
    ylim(0, 3000)
  
  LDSR_POST_DS <- ggpubr::ggscatter(get(paste0('merge_lvl', LEVEL)), x = "LDSR", y = "cnts_ds",
                             add = "reg.line", conf.int = TRUE, 
                             cor.coef = TRUE, cor.method = "pearson",
                             xlab = "LDSR", ylab = "Cell counts", 
                             title = paste0('Cluster level ', LEVEL, ' downsampled')) +
    ylim(0, 3000)
  
  MAGMA_SCATTER <- plot_grid(MAGMA_PRE_DS, MAGMA_POST_DS)
  LDSR_SCATTER <- plot_grid(LDSR_PRE_DS, LDSR_POST_DS)
  
  assign(paste0('magma_dnwSmple_scatter_lvl', LEVEL, '_plot'), MAGMA_SCATTER, envir = .GlobalEnv)  
  assign(paste0('ldsr_dnwSmple_scatter_lvl', LEVEL, '_plot'), LDSR_SCATTER, envir = .GlobalEnv)


}

# Downsampled  
magma_dwnSmple_lvl1_plot <- plot_grid(magma_SCZ_lvl_1_35UP_10DOWN_plot, magma_dwnSmpl_SCZ_lvl_1_35UP_10DOWN_plot)
magma_dwnSmple_lvl2_plot <- plot_grid(magma_SCZ_lvl_2_35UP_10DOWN_plot, magma_dwnSmpl_SCZ_lvl_2_35UP_10DOWN_plot)

# magma_dwnSmple_lvl1_all_plot <- plot_grid(magma_dwnSmpl_SCZ_lvl_1_10UP_10DOWN_plot, magma_dwnSmpl_SCZ_lvl_1_35UP_10DOWN_plot, 
#           magma_dwnSmpl_SCZ_lvl_1_100UP_100DOWN_plot, magma_dwnSmpl_HEIGHT_lvl_1_10UP_10DOWN_plot, 
#           magma_dwnSmpl_HEIGHT_lvl_1_35UP_10DOWN_plot, magma_dwnSmpl_HEIGHT_lvl_1_100UP_100DOWN_plot)

# Gene_windows
magma_windows_lvl1_plot <- plot_grid(magma_SCZ_lvl_1_10UP_10DOWN_plot, magma_SCZ_lvl_1_35UP_10DOWN_plot, 
                                     magma_SCZ_lvl_1_100UP_100DOWN_plot, magma_HEIGHT_lvl_1_10UP_10DOWN_plot, 
                                     magma_HEIGHT_lvl_1_35UP_10DOWN_plot, magma_HEIGHT_lvl_1_100UP_100DOWN_plot, ncol = 3)

ldsr_windows_lvl1_plot <- plot_grid(ldsr_SCZ_lvl_1_10UP_10DOWN_plot, ldsr_SCZ_lvl_1_35UP_10DOWN_plot, 
                                    ldsr_SCZ_lvl_1_100UP_100DOWN_plot, ldsr_HEIGHT_lvl_1_10UP_10DOWN_plot, 
                                    ldsr_HEIGHT_lvl_1_35UP_10DOWN_plot, ldsr_HEIGHT_lvl_1_100UP_100DOWN_plot, ncol = 3)

##  scDRS  ----------------------------------------------------------------------------
scDRS_all_umap <- plot_grid(scDRS_all_SCZ_group_umap, scDRS_all_HEIGHT_group_umap, ncol = 1)
scDRS_MGE_umap <- plot_grid(scDRS_MGE_SCZ_group_umap, scDRS_MGE_HEIGHT_group_umap, ncol = 1)
scDRS_LGE_umap <- plot_grid(scDRS_LGE_SCZ_group_umap, scDRS_LGE_HEIGHT_group_umap, ncol = 1)
scDRS_CGE_umap <- plot_grid(scDRS_CGE_SCZ_group_umap, scDRS_CGE_HEIGHT_group_umap, ncol = 1)
scDRS_Progenitor_umap <- plot_grid(scDRS_Progenitor_SCZ_group_umap, scDRS_Progenitor_HEIGHT_group_umap, ncol = 1)

# Render markdown
rmarkdown::render(MARKDOWN_FILE, output_file = REPORT_FILE, output_dir = REPORT_DIR)

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
