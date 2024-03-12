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
require(readxl)

## Set variables  ---------------------------------------------------------------------
results_dir <- '~/Desktop/fetal_brain_snRNAseq_GE_270922/results/'
resources_dir <- '~/Desktop/fetal_brain_snRNAseq_GE_270922/resources/'
data_dir <- paste0(results_dir, '01R_objects/')
out_dir <- '~/Desktop/fetal_brain_snRNAseq_GE_270922/results/'
genelist_dir <- paste0(out_dir, '02GENE_LISTS/')
scz_dir <- paste0(resources_dir, '')

# monocle_dir <- paste0(genelist_dir, 'diff_exp/')
# diffexp_dir <- paste0(genelist_dir, 'shi_bc/MAGMA_DIFFEXP/')
# ldsr_diffexp_dir <- paste0(genelist_dir, 'shi_bc/LDSR_DIFFEXP/')
# dir.create(paste0(diffexp_dir))
# dir.create(paste0(ldsr_diffexp_dir))
# upstream <- 100000
# downstream <- 100000
# window <- "100UP_100DOWN"

## Load Data  -------------------------------------------------------------------------
seurat_shi <- readRDS(paste0(data_dir, 'seurat_shi_bc.rds'))
prioritsed_genes <- read_excel(paste0(scz_dir, 'Supplementary Table 12.xlsx'), sheet = 'ST12 all criteria') %>%
  select(Symbol.ID, Prioritised) %>%
  filter(Prioritised == 1) %>%
  pull(Symbol.ID)
schema_genes <- read_excel(paste0(scz_dir, 'Supplementary Table 20.xlsx'), sheet = 'Gene Lists') %>%
  select(`SCHEMA genes (FDR < 5%)`) %>%
  filter(!is.na(.)) %>%
  pull()
  

# remove MG
seurat_subset <- subset(x = seurat_shi, subset = cluster_level_1 != "Microglia")

# Convert Seurat obj to SCE
shi_sce <- as.SingleCellExperiment(seurat_subset)

# Run Slingshot
shi_sce <- slingshot(data = shi_sce, clusterLabels = '', reducedDim = 'UMAP', 
                     start.clus = 'Progenitor', end.clus = c('MGE', 'CGE', 'LGE'), 
                     extend = 'n', stretch = 0)

# Plot lineages
plot(reducedDims(shi_sce)$UMAP, col = brewer.pal(9, 'Set1')[shi_sce$cluster_level_1], pch = 16, asp = 1)
lines(SlingshotDataSet(shi_sce), lwd = 2, type = 'lineages', col = 'black')
title('Lineages')


# Plot curves
plot(reducedDims(shi_sce)$UMAP, col = brewer.pal(9,'Set1')[shi_sce$cluster_level_1], pch=16, asp = 1)
lines(SlingshotDataSet(shi_sce), lwd =2, type = 'curves', col = 'black')
title('Curves')

# Check Lineages 
slingLineages(shi_sce)

# L1 - pseudo
colors <- grDevices::colorRampPalette(rev(brewer.pal(11,'Spectral')[-6]))(100)
plotcol <- colors[cut(shi_sce$slingPseudotime_1, breaks=100)]
plot(reducedDims(shi_sce)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(shi_sce), lwd=2, col='black')
title('L1_pseudo')

# L2 - pseudo
plotcol <- colors[cut(shi_sce$slingPseudotime_2, breaks=100)]
plot(reducedDims(shi_sce)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(shi_sce), lwd = 2, col = 'black')
title('L2_pseudo')

# L3 - pseudo
plotcol <- colors[cut(shi_sce$slingPseudotime_3, breaks=100)]
plot(reducedDims(shi_sce)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(shi_sce), lwd=2, col='black')
title('L3_pseudo')

# Plots
tibble(cluster_level_1 = shi_sce@colData$seurat_clusters,
       mge_pseudo = shi_sce@colData$slingPseudotime_1,
       lge_pseudo = shi_sce@colData$slingPseudotime_2,
       cge_pseudo = shi_sce@colData$slingPseudotime_3) %>%
  mutate(all_pseudo = coalesce(mge_pseudo, lge_pseudo,cge_pseudo)) %>%
  ggplot(aes(x = all_pseudo, y = cluster_level_1,
             colour = cluster_level_1)) +
  geom_jitter() + theme_classic() +
  xlab("Slingshot pseudotime") + ylab("cell type") +
  ggtitle("Cells ordered by Slingshot pseudotime")

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
saveRDS(shi_sce_ts, paste0(data_dir, 'sce_shi_ts.rds'))
shi_sce_ts <- readRDS(paste0(data_dir, 'sce_shi.rds'))

# Check for gene convergenece
table(rowData(shi_sce)$tradeSeq$converged)

# Run association test
assoRes <- associationTest(shi_sce)
head(assoRes)

# Run startVsEndTest
# Wald statistic orders results by the degree of difference but no the direction
# LogFC
startRes <- startVsEndTest(shi_sce_ts, lineages = TRUE, l2fc = 0.5)

startRes %>% 
  as_tibble(rownames = 'genes') %>%
  filter(genes %in% c('NEGR1', 'ZNF804A', 'EMX1', 'LINC01088')) %>%
  rename_with(~ gsub("lineage1", "mge", .x, fixed = TRUE)) %>%
  rename_with(~ gsub("lineage2", "lge", .x, fixed = TRUE)) %>%
  rename_with(~ gsub("lineage3", "cge", .x, fixed = TRUE)) %>%
  select(-contains("df")) 

mge_wald_up <- startRes %>%
  as_tibble(rownames = 'genes') %>%
  dplyr::select(genes, waldStat_lineage1, df_lineage1, pvalue_lineage1, logFClineage1) %>%
  filter(logFClineage1 < 0) %>%
  arrange(pvalue_lineage1) %>%
  filter(pvalue_lineage1 < 0.05) 

lge_wald_up <- startRes %>%
  as_tibble(rownames = 'genes') %>%
  dplyr::select(genes, waldStat_lineage2, df_lineage2, pvalue_lineage2, logFClineage2) %>%
  filter(logFClineage2 < 0) %>%
  arrange(desc(logFClineage2)) %>%
  filter(pvalue_lineage2 < 0.05) 

cge_wald_up <- startRes %>%
  as_tibble(rownames = 'genes') %>%
  dplyr::select(genes, waldStat_lineage3, df_lineage3, pvalue_lineage3, logFClineage3) %>%
  filter(logFClineage3 < 0) %>%
  arrange(pvalue_lineage3) %>%
  filter(pvalue_lineage3 < 0.05) 
  
mge_wald_up %>% filter(genes %in% c("NXPH1", "FAM83D", "KLF6", "GPR98"))
lge_wald_up %>% filter(genes %in% c("NXPH1", "FAM83D", "KLF6", "GPR98"))
cge_wald_up %>% filter(genes %in% c("NXPH1", "FAM83D", "KLF6", "GPR98"))


intersect(mge_wald_up %>% pull(genes), prioritsed_genes %>% pull(gene))
intersect(lge_wald_up %>% pull(genes), prioritsed_genes %>% pull(gene))
intersect(cge_wald_up %>% pull(genes), prioritsed_genes %>% pull(gene))

# intersect(mge_wald_up %>% pull(genes), schema_genes)
# intersect(lge_wald_up %>% pull(genes), schema_genes)
# intersect(cge_wald_up %>% pull(genes), schema_genes)
# 
# intersect(mge_wald_up %>% pull(genes), 'NRXN1')
# intersect(lge_wald_up %>% pull(genes), 'NRXN1')
# intersect(cge_wald_up %>% pull(genes), 'NRXN1')

# ge_intersect <- VennDiagram::get.venn.partitions(list(mge = mge_wald_up %>% pull(genes),
#                                                       lge = lge_wald_up %>% pull(genes),
#                                                       cge = cge_wald_up %>% pull(genes)))
# 
# grid::grid.newpage()
# grid::grid.draw(VennDiagram::venn.diagram(list(mge = mge_wald_up %>% pull(genes), 
#                                                lge = lge_wald_up %>% pull(genes),
#                                                cge = cge_wald_up %>% pull(genes)), 
#                                           NULL))
#ge_intersect$..values..$`6`

# Issues - top 2000 genes
# A: Many of the genes with high log2FC changes only have pseudotime values in very few cells
# and often in ribosomal genes or from genes of different lineages (i.e. myeloid) 

plotSmoothers(shi_sce, counts, gene = mge_wald)

oStart <- base::order(startRes$logFClineage1)
sigGeneStart <- names(shi_sce)[oStart[1]]
plotGeneCount(crv, counts(shi_sce_ts), gene = sigGeneStart)

plotSmoothers(shi_sce_ts, counts(shi_sce_ts), gene = 'NPY') + ggtitle('NPY')

plot_list <- list()
for (gene in c("FAM83D", "EMX1", "GPR98", "ZNF804A", "NEGR1", "ZNF804A")) {
  
  plot <- plotSmoothers(shi_sce_ts, counts(shi_sce_ts), gene = gene) + ggtitle(gene)
  plot_list[[gene]] <- plot
  
}
plot_grid(plotlist = plot_list) + ggtitle('test')

plot_list <- list()
for (gene in c('NXPH1', 'FAM83D', 'KLF6', 'GPR98')) {
  
  plot <- plotSmoothers(shi_sce_ts, counts(shi_sce_ts), gene = gene) + ggtitle(gene)
  plot_list[[gene]] <- plot
  
}
plot_grid(plotlist = plot_list, ncol = 2) + ggtitle('test')

plotSmoothers(shi_sce_ts, counts(shi_sce_ts), gene = 'GAD2') + ggtitle('GAD2')
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
