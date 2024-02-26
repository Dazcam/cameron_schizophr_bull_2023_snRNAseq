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

# remove MG
seurat_subset <- subset(x = seurat_shi, subset = cluster_level_1 != "Microglia")

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
saveRDS(shi_sce_ts, paste0(data_dir, 'sce_shi_ts.rds'))

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
plotSmoothers(shi_sce_ts, counts(shi_sce_ts), gene = "FOXP1")
plotGeneCount(shi_sce_ts$crv$pseudotime.Lineage1, counts(shi_sce_ts), gene = sigGeneStart)
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
