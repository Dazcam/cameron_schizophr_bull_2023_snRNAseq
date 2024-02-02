#--------------------------------------------------------------------------------------
#
#    snRNAseq monocle analysis - laptop
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Running monocle analysis using MK167 as the root gene 
#  Using the seurat wrapper code verbatim:
#  http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/monocle3.html

# Load libraries  ---------------------------------------------------------------------
require(monocle3)
require(SeuratWrappers)
require(Seurat)
require(tidyverse)
require(cowplot)

## Set variables  ---------------------------------------------------------------------
results_dir <- '~/Desktop/fetal_brain_snRNAseq_GE_270922/results/'
data_dir <- paste0(results_dir, '01R_objects/')
out_dir <- '~/Desktop/fetal_brain_snRNAseq_GE_270922/results/'
genelist_dir <- paste0(out_dir, '02GENE_LISTS/')
monocle_dir <- paste0(genelist_dir, 'diff_exp/')
diffexp_dir <- paste0(genelist_dir, 'shi_bc/MAGMA_DIFFEXP/')
ldsr_diffexp_dir <- paste0(genelist_dir, 'shi_bc/LDSR_DIFFEXP/')
dir.create(paste0(diffexp_dir))
dir.create(paste0(ldsr_diff_dir))
upstream <- 100000
downstream <- 100000
window <- "100UP_100DOWN"

## Load Data  -------------------------------------------------------------------------
seurat_shi <- readRDS(paste0(data_dir, 'seurat_shi_bc.rds'))

## Load Data  -------------------------------------------------------------------------
cds <- SeuratWrappers::as.cell_data_set(seurat_shi)
cds <- cluster_cells(cds)

## Retrieve info from Monocle class_cell_data object
colData(cds) # Seurat obj meta data
fData(cds)$gene_names <- rownames(fData(cds))
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
cds <- learn_graph(cds, use_partition = FALSE)
cluster_tree_plot <- plot_cells(cds, 
           label_groups_by_cluster = TRUE, 
           label_leaves = TRUE, 
           label_branch_points = TRUE,
           label_roots = TRUE)

# Manually set root based on graph 
cds <- order_cells(cds)

# Visualise
pseudo_plot <- plot_cells(cds, 
                          color_cells_by = 'pseudotime',
                          label_groups_by_cluster = FALSE, 
                          label_leaves = FALSE, 
                          label_branch_points = FALSE,
                          label_roots = TRUE)

orig_cluster_plot <- DimPlot(seurat_shi, group.by = 'seurat_clusters', label = TRUE)
pcw_plot <- DimPlot(seurat_shi, group.by = 'pcw')
merged_cluster_plot <- DimPlot(seurat_shi, group.by = 'cluster_level_1', label = TRUE)

# Visualise pseudotime by cluster - boxplot
cds$pseudotime <- pseudotime(cds)
monocle_meta_tbl <- as_tibble(colData(cds))
pseudo_boxplot <- ggplot(monocle_meta_tbl, aes(pseudotime, reorder(cluster_level_1, pseudotime, median),
                            fill = cluster_level_1)) +
  geom_boxplot()

plot_grid(orig_cluster_plot, merged_cluster_plot, cluster_tree_plot, 
          pseudo_plot, pcw_plot, pseudo_boxplot, ncol = 2)


## Run differential expression analyses
# only.pos = features that are more highly expressed in the ident.1 group.
seurat_early_InN_vs_Prog <- FindMarkers(seurat_shi, ident.1 = "Early_InN", 
                                        ident.2 = "Progenitor", only.pos = TRUE, 
                                        logfc.threshold = 0.58)
seurat_MGE_vs_early_InN <- FindMarkers(seurat_shi, ident.1 = "MGE", 
                                       ident.2 = "Early_InN", only.pos = TRUE, 
                                       logfc.threshold = 0.58)
seurat_CGE_vs_early_InN <- FindMarkers(seurat_shi, ident.1 = "CGE", 
                                       ident.2 = "Early_InN", only.pos = TRUE, 
                                       logfc.threshold = 0.58)
seurat_CGEandLGE_vs_early_InN <- FindMarkers(seurat_shi, ident.1 = c("CGE", "LGE"), 
                                             ident.2 = "Early_InN", only.pos = TRUE, 
                                             logfc.threshold = 0.58)
seurat_LGE_vs_CGE <- FindMarkers(seurat_shi, ident.1 = "LGE", 
                                 ident.2 = "CGE", only.pos = TRUE, 
                                 logfc.threshold = 0.58)

nrow(seurat_early_InN_vs_Prog)
nrow(seurat_CGE_vs_early_InN)
nrow(seurat_MGE_vs_early_InN)
nrow(seurat_CGEandLGE_vs_early_InN)
nrow(seurat_LGE_vs_CGE)

magma_diffExp_list <- list()

# Create gene sets for gene set enrichment analyses
for (conditional in c('early_InN_vs_Prog', 'CGE_vs_early_InN', 
                      'MGE_vs_early_InN','CGEandLGE_vs_early_InN', 
                      'LGE_vs_CGE')) {
  
  entrez_genes <- get(paste0('seurat_', conditional)) %>%
    as_tibble(rownames = 'genes') %>%
    inner_join(gene_coords, by = join_by(genes == hgnc)) %>%
    unique() %>%
    drop_na() %>%
    slice_head(n = 735) %>%
    dplyr::select(entrez) %>%
    with(., split(entrez, conditional))
    #write_tsv(paste0(genelist_dir, 'diff_expr/', conditional, '_diff_exp.tsv'))
    
    magma_diffExp_list <- c(magma_diffExp_list, entrez_genes)
    
}

# Write gene lists to file
for(i in names(magma_diffExp_list)) {
  
  cat(i, " ", paste(magma_diffExp_list[[i]], collapse = " "), "\n", 
      file = paste0(diffexp_dir, 'snRNAseq_GE_diffexp_gene_sets.txt')
      , sep = '', append = TRUE)
  
}
     
for (LINE in seq(1, 4, 1)) {
  
  cond_genelists <- readLines(paste0(diffexp_dir, 'snRNAseq_GE_diffexp_gene_sets.txt'))
  genelist_name <- unlist(strsplit(cond_genelists[as.integer(LINE)], " "))[1]
  
  cat('Obtaining gene coords for', genelist_name, '...\n\n')
  
  df <- unlist(strsplit(cond_genelists[as.integer(LINE)], " ")) %>% 
    as_tibble() %>%
    janitor::row_to_names(row_number = 1) %>%
    dplyr::rename(entrez = 1) %>%
    inner_join(gene_coords) %>%
    dplyr::mutate(start = ifelse(start - upstream < 0, 0, start - upstream), end = end + downstream) %>%
    dplyr::select(chr, start, end, entrez) %>%
    write_tsv(paste0(ldsr_diffexp_dir, genelist_name, '.', window, '.bed'), col_names = FALSE)
  
}

# Find genes that are diff expr across trajectory
deg <- graph_test(cds, neighbor_graph = 'principal_graph', cores = 4)

## Also tried to order cells by whole progenitor population: to manny cells
## And setting cell with highest MKI67 expression to root: better but not best pseudotime
# Use Progenitor population as root
# cds <- order_cells(cds, root_cells = colnames(cds[, clusters(cds) == 'Progenitor']))
# Use MK167 as root
# max.avp <- which.max(unlist(FetchData(seurat_shi, "MKI67")))
# max.avp <- colnames(seurat_shi)[max.avp]
# cds_mk167 <- order_cells(cds, root_cells = max.avp)
# b <- plot_cells(cds_mk167, color_cells_by = "pseudotime", ,
#                 label_groups_by_cluster = TRUE, 
#                 label_leaves = TRUE, 
#                 label_branch_points = TRUE,
#                 label_roots = TRUE)




#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
