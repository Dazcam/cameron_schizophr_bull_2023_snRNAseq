#--------------------------------------------------------------------------------------
#
#    snRNAseq trajectory analysis 
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/monocle3.html
#  https://bioconductor.org/packages/devel/bioc/vignettes/traviz/inst/doc/slingshot.html

### Still to fix fig names

# Load libraries  ---------------------------------------------------------------------
require(monocle3)
require(SeuratWrappers)
require(Seurat)
require(tidyverse)
require(cowplot)
require(tradeSeq)
require(slingshot)
require(RColorBrewer)
require(readxl)

## Set variables  ---------------------------------------------------------------------
results_dir <- '~/Desktop/fetal_brain_snRNAseq_GE_270922/results/'
resources_dir <- '~/Desktop/fetal_brain_snRNAseq_GE_270922/resources/'
data_dir <- paste0(results_dir, '01R_objects/')
scz_dir <- paste0(resources_dir, 'public_data/trubestskoy_2022/2020-08-14908C-s11/')

## Load Data  -------------------------------------------------------------------------
seurat_shi <- readRDS(paste0(data_dir, 'seurat_shi_bc.rds'))
seurat_shi <- subset(x = seurat_shi, subset = cluster_level_1 != "Microglia") # Rm MG

prioritised_genes <- read_excel(paste0(scz_dir, 'Supplementary Table 12.xlsx'), sheet = 'ST12 all criteria') %>%
  select(Symbol.ID, Prioritised) %>%
  filter(Prioritised == 1) %>%
  select(Symbol.ID) %>%
  dplyr::rename(gene = Symbol.ID)

## Load Data  -------------------------------------------------------------------------
cds <- SeuratWrappers::as.cell_data_set(seurat_shi)
cds <- estimate_size_factors(cds) #  https://github.com/cole-trapnell-lab/monocle3/issues/602
cds <- preprocess_cds(cds)  # https://github.com/cole-trapnell-lab/monocle3/issues/655
cds <- cluster_cells(cds)

## Retrieve info from Monocle class_cell_data object
colData(cds) # Seurat obj meta data
fData(cds)$gene_short_name <- rownames(fData(cds))
fData(cds) # Gene names
counts(cds) # raw counts

# Transfer cluster info from seurat object to monocle object
# Set all cells to a single partition
partition_vec <- rep(1, length(cds@colData@rownames))
names(partition_vec) <- cds@colData@rownames
partition_vec <- as.factor(partition_vec)
cds@clusters$UMAP$partitions <- partition_vec

# Transfer cluster info
cds@clusters$UMAP$clusters <- seurat_shi@active.ident

## Run pseudotime analysis  -----------------------------------------------------------
cds <- learn_graph(cds, use_partition = FALSE, close_loop = FALSE)
cluster_tree_plot <- plot_cells(cds, 
                                label_groups_by_cluster = TRUE, 
                                label_leaves = TRUE, 
                                label_branch_points = TRUE,
                                label_roots = TRUE, 
                                label_principal_points = TRUE)

# Manually set root based on graph 
cds <- order_cells(cds, root_pr_nodes = c('Y_125', 'Y_66'))

## Run the graph test -----------------------------------------------------------------
# Default
# Moran's Test - Moran's I give each gene's effect size between -1 and 1
pr_graph_test_res <- graph_test(cds, neighbor_graph = "principal_graph", cores = 3)
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))

# Subset by scz prioritised genes
pr_graph_test_res_tbl <- as_tibble(pr_graph_test_res, rownames = 'gene') %>%
  inner_join(prioritised_genes) %>%
  filter(q_value < 0.05) %>%
  arrange(desc(morans_I))

## Plot results. ----------------------------------------------------------------------
## Fig_S1
cluster_level_1_recode <- colData(cds) %>% 
  as_tibble() %>%
  dplyr::mutate(across(cluster_level_1, str_replace, 'Early_InN', 'IPC')) %>%
  dplyr::mutate(across(cluster_level_1, str_replace, 'GE', 'GE-N')) %>%
  pull(cluster_level_1)
 
cds$cluster_level_1_recode <- cluster_level_1_recode
cds$pseudotime <- pseudotime(cds)

cluster_plot <- plot_cells(cds,
                           color_cells_by = 'cluster_level_1_recode',
                           label_groups_by_cluster = FALSE,
                           label_leaves = FALSE, 
                           label_branch_points = FALSE,
                           label_roots = FALSE, 
                           group_label_size = 4) +
  theme_bw() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        plot.title = element_text(hjust = 0.5, face = 'bold'),
        axis.title.x = element_text(colour = "#000000", size = 14),
        axis.title.y = element_text(colour = "#000000", size = 14),
        axis.text.x  = element_text(colour = "#000000", size = 13, vjust = 0.5),
        axis.text.y  = element_text(colour = "#000000", size = 13)) +
  scale_color_manual(values = c("LGE-N" = "#E69F00", 
                                "MGE-N" = "#56B4E9", 
                                "CGE-N" = "#009E73", 
                                "IPC" = "#F0E442", 
                                "Microglia" = "#0072B2", 
                                "Progenitor" = "#D55E00")) +
  NoLegend()

pseudo_plot <- plot_cells(cds,
                          color_cells_by = 'pseudotime',
                          label_groups_by_cluster = FALSE,
                          label_leaves = FALSE, 
                          label_branch_points = FALSE,
                          label_roots = FALSE) +
  theme_bw() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        plot.title = element_text(hjust = 0.5, face = 'bold'),
        axis.title.x = element_text(colour = "#000000", size = 14),
        axis.title.y = element_text(colour = "#000000", size = 14),
        axis.text.x  = element_text(colour = "#000000", size = 13, vjust = 0.5),
        axis.text.y  = element_text(colour = "#000000", size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_blank()) 
  
# Visualise pseudotime by cluster - boxplot
pseudo_boxplot <- as_tibble(colData(cds)) %>%
  ggplot(aes(pseudotime, 
             reorder(cluster_level_1_recode, pseudotime, median),
             fill = cluster_level_1_recode)) +
  geom_boxplot() +
  scale_fill_manual(values = c("LGE-N" = "#E69F00", 
                               "MGE-N" = "#56B4E9", 
                               "CGE-N" = "#009E73",
                                "IPC" = "#F0E442", 
                                "Progenitor" = "#D55E00")) +
  theme_bw() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        plot.title = element_text(hjust = 0.5, face = 'bold'),
        axis.title.x = element_text(colour = "#000000", size = 14),
        axis.title.y = element_text(colour = "#000000", size = 14),
        axis.text.x  = element_text(colour = "#000000", size = 13, vjust = 0.5),
        axis.text.y  = element_text(colour = "#000000", size = 13)) +
  NoLegend() +
  xlab('Pseudotime') + ylab('Cell Type')

legend <- get_legend(pseudo_plot)
Fig_S1 <- plot_grid(cluster_plot, pseudo_plot + NoLegend(),
                           pseudo_boxplot, ncol = 2, labels = 'AUTO',
                           label_size = 20)
Fig_S1 <- ggdraw(supp_figure_1) +
  draw_plot(legend, .38, .45, .5, .42)


## Fig S3 ----
scz_genes <- c("BCL11B")
lineage_cds <- cds[rowData(cds)$gene_short_name %in% scz_genes,
                   colData(cds)$cluster_level_1_recode %in% c('MGE-N', 'CGE-N', 'LGE-N',
                                                              'Progenitor','IPC')]

lineage_cds$cluster_level_1_recode <- factor(lineage_cds$cluster_level_1_recode, 
                                             levels = c('LGE-N', 'MGE-N', 'CGE-N',
                                                        'IPC', 'Progenitor'))
Fig_S3A <- plot_genes_in_pseudotime(lineage_cds,  
                                    color_cells_by = 'cluster_level_1_recode',
                                    min_expr = 0.5,
                                    horizontal_jitter = T) +
  theme_bw() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        plot.title = element_text(hjust = 0.5, face = 'bold'),
        strip.background = element_rect(fill = 'white', colour = 'white'),
        strip.text.x = element_blank(),
        axis.title.x = element_text(colour = "#000000", size = 14),
        axis.title.y = element_text(colour = "#000000", size = 14),
        axis.text.x  = element_text(colour = "#000000", size = 13, vjust = 0.5),
        axis.text.y  = element_text(colour = "#000000", size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_color_manual(values = c("LGE-N" = "#E69F00", 
                                "MGE-N" = "#56B4E9", 
                                "CGE-N" = "#009E73", 
                                "IPC" = "#F0E442", 
                                "Progenitor" = "#D55E00")) +
 xlab('Pseudotime') + ggtitle(NULL)


## Fig S3B ----
scz_genes <- c("NXPH1")
lineage_cds <- cds[rowData(cds)$gene_short_name %in% scz_genes,
                   colData(cds)$cluster_level_1 %in% c('MGE', 'Progenitor', 'Early_InN')]
lineage_cds$cluster_level_1_recode <- factor(lineage_cds$cluster_level_1_recode, 
                                             levels = c('MGE-N', 'IPC', 'Progenitor'))


Fig_S3B <- plot_genes_in_pseudotime(lineage_cds,  
                                       color_cells_by = 'cluster_level_1_recode',
                                       min_expr = 0.5,
                                       horizontal_jitter = T) +
  theme_bw() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        plot.title = element_text(hjust = 0.5, face = 'bold'),
        strip.background = element_rect(fill = 'white', colour = 'white'),
        strip.text.x = element_blank(),
        axis.title.x = element_text(colour = "#000000", size = 14),
        axis.title.y = element_text(colour = "#000000", size = 14),
        axis.text.x  = element_text(colour = "#000000", size = 13, vjust = 0.5),
        axis.text.y  = element_text(colour = "#000000", size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_color_manual(values = c("LGE-N" = "#E69F00", 
                                "MGE-N" = "#56B4E9", 
                                "CGE-N" = "#009E73", 
                                "IPC" = "#F0E442", 
                                "Microglia" = "#0072B2", 
                                "Progenitor" = "#D55E00")) +
  xlab('Pseudotime') + ggtitle(NULL)

Idents(seurat_shi) <- cluster_level_1_recode
nxp_feature_plot <- FeaturePlot(seurat_shi, features = 'NXPH1', label = T) +
  theme_bw() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        plot.title = element_text(hjust = 0.5, face = 'bold'),
        axis.title.x = element_text(colour = "#000000", size = 14),
        axis.title.y = element_text(colour = "#000000", size = 14),
        axis.text.x  = element_text(colour = "#000000", size = 13, vjust = 0.5),
        axis.text.y  = element_text(colour = "#000000", size = 13),        
        legend.text = element_text(size = 13),
        legend.title = element_blank()) + ggtitle(NULL)

bcl_feature_plot <- FeaturePlot(seurat_shi, features = 'BCL11B', label = T) +
  theme_bw() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 1),
        plot.title = element_text(hjust = 0.5, face = 'bold'),
        axis.title.x = element_text(colour = "#000000", size = 14),
        axis.title.y = element_text(colour = "#000000", size = 14),
        axis.text.x  = element_text(colour = "#000000", size = 13, vjust = 0.5),
        axis.text.y  = element_text(colour = "#000000", size = 13),        
        legend.text = element_text(size = 13),
        legend.title = element_blank()) + ggtitle(NULL)

plot_a <- plot_grid(bcl_curve_plot, bcl_feature_plot)
plot_b <- plot_grid(nxp_curve_plot, nxp_feature_plot)
plot_grid(plot_a, plot_b, ncol = 1, labels = 'AUTO', label_size = 20, 
          align = 'hv', axis = 'tblr')

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------