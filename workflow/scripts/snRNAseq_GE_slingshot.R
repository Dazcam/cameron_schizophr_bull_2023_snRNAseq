#--------------------------------------------------------------------------------------
#
#    snRNAseq trajectory analysis - Hawk
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Run trajectory inference analyses using slingshot / tradeseq on Hawk
#  K and fitGAM functions are resource intensive so I can't run it locally
#  This is run in container

# Load libraries  ---------------------------------------------------------------------
message('Loading Libraries and setting env variables...')
require(Seurat)
require(tradeSeq)
require(slingshot)
require(BiocParallel)

## Set variables  ---------------------------------------------------------------------
results_dir <- toString(snakemake@params[['results_dir']])
data_dir <- paste0(results_dir, '01R_objects/')

# Set resources
message('Setting seed and workers ...')
set.seed(10) # fitGAM is stochastic
message('Cores avaialable on Hawk: ', parallel::detectCores())
BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers <- 20 # use 3 cores
BPPARAM # lists current options

## Load Data  -------------------------------------------------------------------------
message('Loading seurat data ...')
seurat_shi <- readRDS(paste0(data_dir, 'seurat_shi_bc.rds'))

# remove MG
message('Removing MG cells ...')
seurat_subset <- subset(x = seurat_shi, subset = cluster_level_1 != "Microglia")

# Convert Seurat obj to SCE
message('Converting seurat obj to SCE obj ...')
shi_sce <- as.SingleCellExperiment(seurat_subset)

# Run Slingshot
message('Running trajectory inference with Slingshot ...')
shi_sce <- slingshot(data = shi_sce, clusterLabels = 'cluster_level_1', reducedDim = 'UMAP', 
                     start.clus = 'Progenitor', end.clus = c('MGE', 'CGE', 'LGE'))

# Check Lineages 
message('Lineages identities are: ...')
slingLineages(shi_sce)

# Get weights and pseudotime
message('Extract cell weights and pseudotime from SCE object ...')
weights <- get_cell_weights(shi_sce)
pseudotime <- get_pseudotime(shi_sce, weights = weights)

# Get Knots
message('Getting knots ...')
aicMat <- evaluateK(counts = counts(shi_sce), pseudotime = pseudotime, cellWeights = cellWeights
                    , k = 3:20, nGenes = 500, verbose = TRUE, plot = TRUE, parallel = TRUE)

# Run fitGAM
message('Running fitGAM ...')
shi_sce <- fitGAM(counts = counts(shi_sce), pseudotime = pseudotime, cellWeights = cellWeights,
                  nknots = 7, verbose = FALSE, parallel = T, genes = var_genes)

# Save RDS
message('Writing SCE object ...')
saveRDS(shi_sce, paste0(data_dir, 'sce_shi.rds'))

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------