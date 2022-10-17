#--------------------------------------------------------------------------------------
#
#    Shi 2021 snRNAseq data - Seurat testing
#
#--------------------------------------------------------------------------------------

##  Info  -----------------------------------------------------------------------------

# Seurat testing for Shi et al (2021) data
# Tested Shi parameters from methods as best I could - all clusters
# Tested subset of data removing ExNs with / without batch correction 
# Testing res levels: 'snRNAseq_GE_functions.R'

##  Load Packages  --------------------------------------------------------------------
if (!require("Require")) install.packages("Require")
Require::Require(c("tidyverse", "readxl", "data.table", "ggdendro", "Seurat", 
                   "SeuratWrappers", "cowplot", "scCustomize", "rmarkdown")) # scCustomize for contrast in cols
# BiocManager::install(c("EWCE", "AnnotationDbi", "org.Hs.eg.db", "scuttle", zellkonverter"))

##  Set Variables  --------------------------------------------------------------------
DATA_DIR <- '~/Desktop/fetal_brain_snRNAseq_GE_270922/resources/'
OUT_DIR <- '~/Desktop/fetal_brain_snRNAseq_GE_270922/results/'
SCRIPT_DIR <- '~/Desktop/fetal_brain_snRNAseq_GE_270922/workflow/scripts/'
MARKDOWN_FILE <- paste0(SCRIPT_DIR, 'snRNAseq_GE_seurat.Rmd')
REPORT_DIR <- paste0(OUT_DIR, 'rmarkdown_reports/')
REPORT_FILE <- 'snRNAseq_GE_seurat_testing.html'

##  Load scripts an functions   -------------------------------------------------------
source(paste0(SCRIPT_DIR, 'snRNAseq_GE_prep_shi_data_for_Seurat.R'))
source(paste0(SCRIPT_DIR, 'snRNAseq_GE_functions.R'))
source(paste0(SCRIPT_DIR, 'snRNAseq_GE_marker_genes.R'))

## Replicate Shi results based on the parameters specified in Methods -----------------
# Uses fastMNN for batch correction - from SeuratWrappers
# Assumption that MNN was set on pcw / batch
seurat.mnn <- CreateSeuratObject(counts = shi_data, meta.data = shi_meta)
seurat.mnn <- subset(seurat.mnn, subset = pcw %in% c("09", "12_01", "13", "16", "18_01", 
                                                           "12_02", "12_02.1"))
seurat.mnn <- NormalizeData(seurat.mnn)
seurat.mnn  <- FindVariableFeatures(seurat.mnn, selection.method = "vst", nfeatures = 2000)
seurat.mnn  <- RunFastMNN(object.list = SplitObject(seurat.mnn, split.by = "pcw"))
seurat.mnn  <- RunUMAP(seurat.mnn , reduction = "mnn", dims = 1:10)
seurat.mnn  <- FindNeighbors(seurat.mnn , reduction = "mnn", dims = 1:10)
seurat.mnn  <- FindClusters(seurat.mnn, resolution = c(0.2, 0.3, 0.4, 0.5))

# Reporting
seurat_resolution_test(seurat.mnn, 
                       c(0.2, 0.3, 0.4, 0.5),
                       'pcw',
                       shi1C_markers,
                       'shi_raw')


## Subset Shi data - no batch correction  ---------------------------------------------
seurat.shi <- CreateSeuratObject(counts = shi_data, meta.data = shi_meta)
seurat.shi <- subset(seurat.shi, subset = pcw %in% c("09", "12_01", "13", "16", "18_01", 
                                                           "12_02", "12_02.1"))
seurat.shi <- subset(seurat.shi, subset = ClusterID %in% c("MGE","LGE", "CGE", "OPC", "Microglia",
                                                           "Endothelial", "progenitor"))
seurat.shi <- NormalizeData(seurat.shi)
seurat.shi <- FindVariableFeatures(seurat.shi, selection.method = "vst", nfeatures = 2000) 
seurat.shi <- ScaleData(seurat.shi) 
seurat.shi <- RunPCA(seurat.shi, features = VariableFeatures(object = seurat.shi)) 
seurat.shi <- FindNeighbors(seurat.shi) # default 
seurat.shi <- FindClusters(seurat.shi, resolution = c(0.2, 0.3, 0.4, 0.5)) # Resolution default is 0.8
seurat.shi <- RunUMAP(seurat.shi, dims = 1:10)

# Reporting
seurat_resolution_test(seurat.shi, 
                       c(0.2, 0.3, 0.4, 0.5),
                       'pcw',
                       our_markers,
                       'shi_subset')


## Subset Shi data - with batch correction  ---------------------------------------------
seurat.shi.bc <- CreateSeuratObject(counts = shi_data, meta.data = shi_meta)
seurat.shi.bc <- subset(seurat.shi.bc, subset = pcw %in% c("09", "12_01", "13", "16", "18_01", 
                                                           "12_02", "12_02.1"))
seurat.shi.bc <- subset(seurat.shi.bc, subset = ClusterID %in% c("MGE","LGE", "CGE", "OPC", "Microglia",
                                                           "Endothelial", "progenitor"))
seurat.shi.bc <- NormalizeData(seurat.shi.bc)
seurat.shi.bc  <- FindVariableFeatures(seurat.shi.bc, selection.method = "vst", nfeatures = 2000)
seurat.shi.bc  <- RunFastMNN(object.list = SplitObject(seurat.shi.bc, split.by = "pcw"))
seurat.shi.bc  <- RunUMAP(seurat.shi.bc , reduction = "mnn", dims = 1:10)
seurat.shi.bc  <- FindNeighbors(seurat.shi.bc , reduction = "mnn", dims = 1:10)
seurat.shi.bc  <- FindClusters(seurat.shi.bc, resolution = c(0.2, 0.3, 0.4, 0.5))

# Reporting
seurat_resolution_test(seurat.shi.bc, 
                       c(0.2, 0.3, 0.4, 0.5),
                       'pcw',
                       our_markers,
                       'shi_subset_bc')

cluster.markers <- FindMarkers(seurat.shi.bc, ident.1 = 3, min.pct = 0.25)

# Create markdown file
render(MARKDOWN_FILE, output_file = REPORT_FILE, output_dir = REPORT_DIR)


## Subset Shi data - remove pcw clusters  ---------------------------------------------
seurat.shi.bc <- CreateSeuratObject(counts = shi_data, meta.data = shi_meta)
seurat.shi.bc <- subset(seurat.shi.bc, subset = pcw %in% c("09", "12_01", "13", "16", "18_01", 
                                                           "12_02", "12_02.1"))
seurat.shi.bc <- subset(seurat.shi.bc, subset = ClusterID %in% c("MGE","LGE", "CGE", "OPC", "Microglia",
                                                                 "Endothelial", "progenitor"))
seurat.shi.bc <- NormalizeData(seurat.shi.bc)
seurat.shi.bc  <- FindVariableFeatures(seurat.shi.bc, selection.method = "vst", nfeatures = 2000)
seurat.shi.bc  <- RunFastMNN(object.list = SplitObject(seurat.shi.bc, split.by = "pcw"))
seurat.shi.bc  <- RunUMAP(seurat.shi.bc , reduction = "mnn", dims = 1:10)
seurat.shi.bc  <- FindNeighbors(seurat.shi.bc , reduction = "mnn", dims = 1:10)
seurat.shi.bc  <- FindClusters(seurat.shi.bc, resolution = c(0.2, 0.3, 0.4, 0.5))

# Plotting
basic_plot <- DimPlot_scCustom(seurat.shi.bc,
                               DiscretePalette_scCustomize(num_colors = 26, 
                                                           palette = "ditto_seq"))
pcw_plot <- DimPlot_scCustom(seurat.shi.bc, group.by = 'pcw')
shi_plot <- DimPlot_scCustom(seurat.shi.bc, group.by = 'ClusterID')
group_plot <- plot_grid(basic_plot, pcw_plot, shi_plot)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------


## Testing code -----------------------------------------------------------------------
human_colors_list <- c("firebrick1", "dodgerblue", "navy", "forestgreen", "darkorange2", 
                       "darkorchid3", "orchid", "orange", "gold", "gray", "black")
Stacked_VlnPlot(seurat_object = seurat.mnn, features = umap_markers, x_lab_rotate = TRUE,
                colors_use = human_colors_list, group.by = "RNA_snn_res.0.3")
DotPlot_scCustom(seurat_object = seurat.mnn, features = umap_markers, 
                 colors_use = viridis_plasma_dark_high, flip_axes = TRUE)


## Testing code
seurat.test <- CreateSeuratObject(counts = shi_data, meta.data = shi_meta)
seurat.test <- NormalizeData(seurat.test)
seurat.test <- subset(seurat.test, subset = ClusterID %in% c("MGE","LGE", "CGE", "OPC", "Microglia",
                                                           "Endothelial", "progenitor"))
seurat.test <- FindVariableFeatures(seurat.test, selection.method = "vst", nfeatures = 2000) 
seurat.test <- ScaleData(seurat.test) 
seurat.test <- RunPCA(seurat.test, features = VariableFeatures(object = seurat.test)) 
##  Check total variance in 1st 50pcs of data  ----------------------------------------
mat <- GetAssayData(seurat.test , assay = "RNA", slot = "scale.data")
pca <- seurat.test[["pca"]]

# Get the total variance (per row in scaled GeX matrix)
total_variance <- sum(matrixStats::rowVars(mat))

eigValues = (pca@stdev) ^ 2  ## EigenValues
varExplained = eigValues / total_variance

# How much of total variance is captured in 1st 50PCs
sum(varExplained)




## QCs
VlnPlot(seurat.test, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(seurat.test, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat.test, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

head(VariableFeatures(seurat.test), 50)
print(seurat.test[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(seurat.test, reduction = "pca", dims = c(4,3))
DimHeatmap(seurat.test, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(seurat.test, ndims = 50)

seurat.test <- FindNeighbors(pbmc, dims = 1:10)
seurat.test <- FindClusters(seurat.test) # Resolution default is 0.8
seurat.test <- RunUMAP(seurat.test, dims = 1:10)





# Calculate DeX genes
Idents(object = seurat.shi) <- seurat.shi@meta.data$'RNA_snn_res.0.3'
diff_eX_genes <- FindAllMarkers(seurat.shi)

# Display GeX - https://github.com/satijalab/seurat/issues/3211
diff_eX_genes %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(seurat.mnn, features = top10$gene, slot = 'data') + NoLegend()