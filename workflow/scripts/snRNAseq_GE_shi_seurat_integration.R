#--------------------------------------------------------------------------------------
#
#    Shi 2021 snRNAseq data - Integration draft
#
#--------------------------------------------------------------------------------------

##  Info  -----------------------------------------------------------------------------

# Integrating Shi et al and Cameron et al data

##  Load Packages  --------------------------------------------------------------------
if (!require("Require")) install.packages("Require")
Require::Require(c("tidyverse", "readxl", "data.table", "ggdendro", "Seurat", 
                   "SeuratWrappers", "cowplot", "scCustomize", "rmarkdown")) 

##  Set Variables  --------------------------------------------------------------------
DATA_DIR <- '~/Desktop/fetal_brain_snRNAseq_110122/resources/' 
SCRIPT_DIR <- '~/Desktop/fetal_brain_snRNAseq_GE_270922/workflow/scripts/'

##  Load scripts and functions   -------------------------------------------------------
source(paste0(SCRIPT_DIR, 'snRNAseq_GE_prep_shi_data_for_Seurat.R'))
source(paste0(SCRIPT_DIR, 'snRNAseq_GE_functions.R'))
source(paste0(SCRIPT_DIR, 'snRNAseq_GE_marker_genes.R'))

## Load and prep data  ----------------------------------------------------------------
seurat.wge <- readRDS(paste0(DATA_DIR, 'R_objects/seurat.wge.final.rds'))
seurat.wge <- DietSeurat(seurat.wge) # Remove scale data / clean object
seurat.wge$pcw <- rep('15', ncol(seurat.wge))
seurat.wge$study <- rep('cameron', ncol(seurat.wge))

seurat.shi <- CreateSeuratObject(counts = shi_data, meta.data = shi_meta)
seurat.shi <- subset(seurat.shi, subset = pcw %in% c("09", "12_01", "13", "16", "18_01", 
                                                           "12_02", "12_02.1"))
seurat.shi <- subset(seurat.shi, subset = ClusterID %in% c("MGE","LGE", "CGE", "OPC", "Microglia",
                                                           "Endothelial", "progenitor"))

seurat.shi$study <- rep('shi', ncol(seurat.shi))


##  Run integration using scTransform  ------------------------------------------------
# https://github.com/satijalab/seurat/issues/1720
obj.list <- list(seurat.shi, seurat.wge)
obj.list <- lapply(X = obj.list, FUN = function(x) {
  x <- SCTransform(x)
  # Don't run ScaleData here
  x <- RunPCA(x, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 2000)
obj.list <- PrepSCTIntegration(object.list = obj.list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features, 
                                  normalization.method = "SCT", reduction = "rpca")
seurat.int <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

# Don't run ScaleData here after integration
seurat.int <- RunPCA(seurat.int, verbose = FALSE)
seurat.int <- RunUMAP(seurat.int, dims = 1:10)
seurat.int <- FindNeighbors(seurat.int, reduction = "pca", dims = 1:10)
seurat.int <- FindClusters(seurat.int, resolution = 0.5)

# Reporting
seurat_resolution_test(seurat.int, 
                       c(0.2, 0.3, 0.4, 0.5),
                       'pcw',
                       our_markers,
                       'shi_int')

# Plotting
basic_plot <- DimPlot_scCustom(seurat.int,
                               DiscretePalette_scCustomize(num_colors = 26, 
                                                           palette = "ditto_seq"))
pcw_plot <- DimPlot_scCustom(seurat.int, group.by = 'pcw')
study_plot <- DimPlot_scCustom(seurat.int, group.by = 'study')
shi_plot <- DimPlot_scCustom(seurat.int, group.by = 'ClusterID')
cam_plot <- DimPlot_scCustom(seurat.int, group.by = 'cellIDs')

group_plot <- plot_grid(basic_plot, pcw_plot, study_plot, shi_plot, cam_plot)

# Batch correction
seurat.int <- RunFastMNN(object.list = SplitObject(seurat.int, split.by = "study"))
# Throws error: issue with SummarizedExperiment package
# See here: https://github.com/satijalab/seurat/issues/5329

# Error in SummarizedExperiment::SummarizedExperiment(assays = assays) : 
#   the rownames and colnames of the supplied assay(s) must be NULL or identical to those of
# the SummarizedExperiment object (or derivative) to construct


#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------