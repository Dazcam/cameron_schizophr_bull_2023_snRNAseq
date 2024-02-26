#--------------------------------------------------------------------------------------
#
#    snRNAseq trajecory analysis 
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/monocle3.html
#  https://bioconductor.org/packages/devel/bioc/vignettes/traviz/inst/doc/slingshot.html

# Load libraries  ---------------------------------------------------------------------
require(monocle3)
require(SeuratWrappers)
require(Seurat)
require(tidyverse)
require(cowplot)
require(tradeSeq)
require(slingshot)
require(RColorBrewer)

## Set variables  ---------------------------------------------------------------------
results_dir <- '~/Desktop/fetal_brain_snRNAseq_GE_270922/results/'
data_dir <- paste0(results_dir, '01R_objects/')
out_dir <- '~/Desktop/fetal_brain_snRNAseq_GE_270922/results/'
genelist_dir <- paste0(out_dir, '02GENE_LISTS/')
monocle_dir <- paste0(genelist_dir, 'diff_exp/')
diffexp_dir <- paste0(genelist_dir, 'shi_bc/MAGMA_DIFFEXP/')
ldsr_diffexp_dir <- paste0(genelist_dir, 'shi_bc/LDSR_DIFFEXP/')
dir.create(paste0(diffexp_dir))
dir.create(paste0(ldsr_diffexp_dir))
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
                                #trajectory_graph_color = "#d0fefe",
                                label_groups_by_cluster = TRUE, 
                                label_leaves = TRUE, 
                                label_branch_points = TRUE,
                                label_roots = TRUE, 
                                label_principal_points = TRUE)

# Manually set root based on graph 
cds <- order_cells(cds, root_pr_nodes = c('Y_21', 'Y_24', 'Y_137'))



## Add Velmashev code  ----------------------------------------------------------------

#run monocle 3
shi_counts_mat <- GetAssayData(object = seurat_shi, assay = "RNA", slot = "counts")
shi_genes_df = as.data.frame(rownames(shi_counts))
rownames(shi_genes_df) = rownames(shi_counts_mat)
colnames(shi_genes_df) = 'gene_short_name'
shi_meta = seurat_shi@meta.data
cds_vel = new_cell_data_set(shi_counts_mat, cell_metadata = shi_meta, gene_metadata = shi_genes_df)
s.umap <- seurat_shi@"reductions"$umap[[]]
s.umap = s.umap[colnames(cds_vel),]
reducedDims(cds_vel)$"UMAP" <- s.umap
s.clusters = as.character(Idents(seurat_shi))
names(s.clusters) <- names(Idents(seurat_shi))
s.clusters = s.clusters[colnames(cds_vel)]
cds_vel@clusters$"UMAP"$"clusters" <- s.clusters
cds_vel@clusters$UMAP$partitions <- cds_vel@clusters$UMAP$clusters
cds_vel@clusters$UMAP$partitions[cds@clusters$UMAP$partitions != "1"] <- "1"
cds_vel <- learn_graph(cds_vel, use_partition = FALSE)
cds_test <- import_monocle(cds_vel) # Adds graphs to plot

###isolate lineages
###isolate lineages
# Progenitors
lineage = "Prog_to_IPC"
start = 20
end = 112
inc.node = c("Y_135")
cds_test <- isolate_graph(cds_test, start, end, lineage, include_nodes = inc.node)
sel.cluster = c("Progenitor", "Early_InN")
cd_test <- isolate_lineage(cds_test, lineage, sel_clusters = sel.cluster, cl = 4, N = 5)




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
     
for (LINE in seq(1, 5, 1)) {
  
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

# Run gene modules
# This needed to be run first: 
# cds <- preprocess_cds(cds) # May need to add this to final script
gene_module_df <- find_gene_modules(cds[diff_genes,], resolution=1e-2)

test <-gene_module_df %>%
#  group_by(module) %>%
#  count() %>%
  filter(module == 62)
  summary()

cell_group_df <- tibble::tibble(cell = row.names(colData(cds)), 
                                cell_group = clusters(cds)[colnames(cds)])
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
#colnames(agg_mat) <- stringr::str_c("Cluster ", colnames(agg_mat))

pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="column", clustering_method="ward.D2",
                   fontsize=6)

plot_cells(cds,
           genes=gene_module_df %>% 
             dplyr::filter(module %in% c(62, 14, 23, 67, 75, 1, 68)),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

CGE_genes <- c("SCGN", "PROX1", "NFIB", "NFIX", "CALB2")
CGE_lineage_cds <- cds[rowData(cds)$gene_short_name %in% CGE_genes,
                       colData(cds)$cluster_level_1 %in% c("CGE")]
CGE_lineage_cds <- order_cells(CGE_lineage_cds, root_pr_nodes = c('Y_21', 'Y_24', 'Y_137'))
plot_genes_in_pseudotime(CGE_lineage_cds,
                         color_cells_by="pseudotime",
                         min_expr=0.5)

cds_root <- get_graph_root(cds, 'cluster_level_1', 'Progenitor')

# Check that
cds_Y55 <- order_cells(cds_Y45, root_pr_nodes = c('Y_55'), )
pseudo_plot <- plot_cells(cds_Y55, 
                          color_cells_by = 'pseudotime',
                          label_groups_by_cluster = FALSE, 
                          label_leaves = FALSE, 
                          label_branch_points = FALSE,
                          label_roots = TRUE,
                          label_principal_points = TRUE)


# Learn graph for a subset of the data
start_node <- 'Y_61'
end_nodes <- c('Y_5', 'Y_119')
cell_select <- choose_graph_segments(cds = cds, starting_pr_node = start_node, ending_pr_nodes = end_nodes, return_list = TRUE)

cds_mge <- cds
cds_mge$cell_subset <- rownames(cds_mge@colData) %in% cell_select$cells 
# Remove MG cells
mge_cell_ids <- cds_mge@colData %>%
  as_tibble(rownames = 'cell_id') %>%
  filter(cluster_level_1 == 'Microglia') %>%
  dplyr::pull(cell_id)

cell_select$cells <- cell_select$cells[ !cell_select$cells %in% mge_cell_ids ]
cds_mge <- cds_mge[,cell_select$cells]
colData(cds_mge)[,cell_select$cells]
cds_mge <- learn_graph(cds_mge, use_partition = FALSE)
order_cells(cds_mge)

# remove MG
seurat_subset <- subset(x = seurat_shi, subset = cluster_level_1 != "Microglia")

# Run Slingshot
# Convert Seurat obj to SCE
shi_sce <- as.SingleCellExperiment(seurat_subset)

# Run Slingshot
shi_sce <- slingshot(data = shi_sce, clusterLabels = 'cluster_level_1', reducedDim = 'UMAP', 
                     start.clus = 'Progenitor', end.clus = c('MGE', 'CGE', 'LGE'))

# Plot lineages
plot(reducedDims(shi_sce)$UMAP, col = brewer.pal(9,'Set1')[shi_sce$cluster_level_1], pch=16, asp = 1)
lines(SlingshotDataSet(shi_sce), lwd = 2, type = 'lineages', col = 'black')
title('Lineages')

# Plot curves
plot(reducedDims(shi_sce)$UMAP, col = brewer.pal(9,'Set1')[shi_sce$cluster_level_1], pch=16, asp = 1)
lines(SlingshotDataSet(shi_sce), lwd =2, type = 'curves', col = 'black')
title('Curves')

# Check Lineages 
slingLineages(shi_sce)

# L1 - pseudo
colors <- grDevices::colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(shi_sce$slingPseudotime_1, breaks=100)]
plot(reducedDims(shi_sce)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(shi_sce), lwd=2, col='black')
title('L1_pseudo')

# L2 - pseudo
plotcol <- colors[cut(shi_sce$slingPseudotime_2, breaks=100)]
plot(reducedDims(shi_sce)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(shi_sce), lwd=2, col='black')
title('L2_pseudo')

# L3 - pseudo
plotcol <- colors[cut(shi_sce$slingPseudotime_3, breaks=100)]
plot(reducedDims(shi_sce)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(shi_sce), lwd=2, col='black')
title('L3_pseudo')


# Run fitGAM
set.seed(6) # fitGAM is stochastic
BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers <- 3 # use 3 cores
BPPARAM # lists current options

# Pull out pseudotime and cellweight values
pseudotime <- slingPseudotime(shi_sce, na = FALSE)
cellWeights <- slingCurveWeights(shi_sce)

# Get K - 7.5 hours to run locally on 3 workers for 5 knots (they recommend 3-8 knots)
aicMat <- evaluateK(counts = counts(shi_sce), pseudotime = pseudotime, cellWeights = cellWeights, 
                    k = 3:8, nGenes = 100, verbose = TRUE, plot = TRUE, parallel = TRUE)
plot_evalutateK_results(aicMat, 3:12)

print(aicMat[1:2, ])


# Can take days to run - use variable genes
# https://github.com/statOmics/tradeSeq/issues/220

# Pull out var genes to improve fitGAM running time
var_genes <- seurat_shi@assays$mnn.reconstructed@var.features

# Run fit GAM - takes 2 days locally (3 workers) using 7 knots and 2000 var genes
# Note that subsetting sce object to top 2000 genes will omit 
shi_sce_ts <- fitGAM(counts = counts(shi_sce), pseudotime = pseudotime, cellWeights = cellWeights,
                  nknots = 7, verbose = FALSE, parallel = T, genes = var_genes)

# Save sce object
saveRDS(shi_sce, paste0(data_dir, 'sce_shi.rds'))

# Check for gene convergenece
table(rowData(shi_sce)$tradeSeq$converged)

# Run association test
assoRes <- associationTest(shi_sce)
head(assoRes)

startRes <- startVsEndTest(shi_sce_ts, lineages = TRUE, l2fc = log2(2))

mge_wald_up <- startRes %>%
  as_tibble(rownames = 'genes') %>%
  dplyr::select(genes, waldStat_lineage1, df_lineage1, pvalue_lineage1, logFClineage1) %>%
  filter(logFClineage1 > 0) %>%
  arrange(desc(logFClineage1))

lge_wald_up <- startRes %>%
  as_tibble(rownames = 'genes') %>%
  dplyr::select(genes, waldStat_lineage2, df_lineage2, pvalue_lineage2, logFClineage2) %>%
  filter(logFClineage2 > 0) %>%
  arrange(desc(logFClineage2))

cge_wald_up <- startRes %>%
  as_tibble(rownames = 'genes') %>%
  dplyr::select(genes, waldStat_lineage3, df_lineage3, pvalue_lineage3, logFClineage3) %>%
  filter(logFClineage3 > 0) %>%
  arrange(desc(logFClineage3))


length(intersect(mge_wald_up %>% pull(genes), lge_wald_up %>% pull(genes)))
length(intersect(mge_wald_up %>% pull(genes), cge_wald_up %>% pull(genes)))
length(intersect(lge_wald_up %>% pull(genes), cge_wald_up %>% pull(genes)))
length(Reduce(intersect, list(lge_wald_up %>% pull(genes), 
                              cge_wald_up %>% pull(genes),
                              mge_wald_up %>% pull(genes))))

ge_intersect <- VennDiagram::get.venn.partitions(list(mge = mge_wald_up %>% pull(genes),
                                                      lge = lge_wald_up %>% pull(genes),
                                                      cge = cge_wald_up %>% pull(genes)))

grid::grid.newpage()
grid::grid.draw(VennDiagram::venn.diagram(list(mge = mge_wald_up %>% pull(genes), 
                                               lge = lge_wald_up %>% pull(genes),
                                               cge = cge_wald_up %>% pull(genes)), 
                                          NULL))


ge_intersect$..values..$`6`

# Issues - top 2000 genes
# A: Many of the genes with high log2FC changes only have pseudotime values in very few cells
# and often in ribosomal genes or from genes of different lineages (i.e. myeloid) 

plotSmoothers(shi_sce, counts, gene = mge_wald)

oStart <- base::order(startRes$logFClineage1)
sigGeneStart <- names(shi_sce_ts)[oStart[1]]
plotSmoothers(shi_sce_ts, counts(shi_sce_ts), gene = "GAPDH")
plotGeneCount(shi_sce_ts$crv$pseudotime.Lineage1, counts(shi_sce_ts), gene = sigGeneStart)
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
