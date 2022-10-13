#--------------------------------------------------------------------------------------
#
#    Shi 2021 snRNAseq data - Seurat testing
#
#--------------------------------------------------------------------------------------

##  Info  --------------------------------------------------------------------

# 1. Seurat testing for Shi et al (2021) data

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

## Markers
mge_markers <- c("PLS3", "NXPH1", "SFTA3", "SOX6")
lge_markers <- c("SIX3", "ZNF503", "SERTAD4", "ISL1")
cge_markers <- c("NFIB", "PDZRN3", "AP1S2", "CALB2", "SCGN",
                 "PCDH9", "KLHL35", "ANKS1B")
supp1C_markers <- c("HMGN2", "EOMES", "CD68", "OLIG1", "OLIG2", 
                    "IGFBP7", "LHX2", "GBX2", "LHX9")
our_markers <- c('GAD1', 'GAD2', 'SLC32A1', 'GLI3', 'TNC',  
                 'SLC17A7', 'SPI1', 'PROX1', 'SCGN', 'LHX6', 
                 'LHX8', 'NXPH1', 'MEIS2', 'ZFHX3', 'ISL1')
shi1C_markers <- c('HES1', 'MKI67', 'ASCL1', 'NKX2-1', 'LHX6',
                   'SOX6', 'NR2F1', 'NR2F2', 'FOXP1', 'MEIS2', 
                   'NEUROD2', 'TCF7L2')


## Replicate Shi results based on the parameters specified in Methods -----------------
# Uses fastMNN for batch correction - from SeuratWrappers
# Assumption that MNN was set on pcw / batch
seurat.mnn <- CreateSeuratObject(counts = shi_data, meta.data = shi_meta)
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
seurat.shi.bc <- subset(seurat.shi.bc, subset = ClusterID %in% c("MGE","LGE", "CGE", "OPC", "Microglia",
                                                           "Endothelial", "progenitor"))
seurat.shi.bc <- NormalizeData(seurat.shi.bc)
seurat.shi.bc  <- FindVariableFeatures(seurat.shi.bc, selection.method = "vst", nfeatures = 2000)
seurat.shi.bc  <- RunFastMNN(object.list = SplitObject(seurat.shi.bc, split.by = "pcw"))
seurat.shi.bc  <- RunUMAP(seurat.shi.bc , reduction = "mnn", dims = 1:10)
seurat.shi.bc  <- FindNeighbors(seurat.shi.bc , reduction = "mnn", dims = 1:10)
seurat.shi.bc  <- FindClusters(seurat.shi.bc, resolution = c(0.2, 0.3, 0.4, 0.5))

seurat_resolution_test(seurat.shi.bc, 
                       c(0.2, 0.3, 0.4, 0.5),
                       'pcw',
                       our_markers,
                       'shi_subset_bc')

# Create markdown file
render(MARKDOWN_FILE, output_file = REPORT_FILE, output_dir = REPORT_DIR)






## Testing code -----------------------------------------------------------------------
human_colors_list <- c("firebrick1", "dodgerblue", "navy", "forestgreen", "darkorange2", "darkorchid3", "orchid",
                       "orange", "gold", "gray", "black")
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