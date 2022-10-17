#--------------------------------------------------------------------------------------
#
#    Shi 2021 snRNAseq data - Seurat 
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

##  Load scripts and functions   -------------------------------------------------------
source(paste0(SCRIPT_DIR, 'snRNAseq_GE_prep_shi_data_for_Seurat.R'))
source(paste0(SCRIPT_DIR, 'snRNAseq_GE_functions.R'))
source(paste0(SCRIPT_DIR, 'snRNAseq_GE_marker_genes.R'))

## Subset Shi data - with batch correction  ---------------------------------------------
seurat.shi <- CreateSeuratObject(counts = shi_data, meta.data = shi_meta)
seurat.shi.bc <- subset(seurat.shi, subset = pcw %in% c("09", "12_01", "13", "16", "18_01", 
                                                           "12_02", "12_02.1"))
seurat.shi.bc <- subset(seurat.shi.bc, subset = ClusterID %in% c("MGE","LGE", "CGE", "OPC", "Microglia",
                                                                 "Endothelial", "progenitor"))
seurat.shi.bc <- NormalizeData(seurat.shi.bc)
seurat.shi.bc <- FindVariableFeatures(seurat.shi.bc, selection.method = "vst", nfeatures = 2000)
seurat.shi.bc <- RunFastMNN(object.list = SplitObject(seurat.shi.bc, split.by = "pcw"))
seurat.shi.bc <- RunUMAP(seurat.shi.bc, reduction = "mnn", dims = 1:10)
seurat.shi.bc <- FindNeighbors(seurat.shi.bc, reduction = "mnn", dims = 1:10)
seurat.shi.bc <- FindClusters(seurat.shi.bc, resolution = 0.4)

# Rename clusters
new_idents <- c('MGE', 'LGE', 'CGE', 'Early InN', 'Progenitor', 
                'Progenitor', 'Progenitor', 'MGE', 'Progenitor', 'Progenitor',
                'Microglia', 'Endothelial')
seurat.shi.bc <- Rename_Clusters(seurat.shi.bc, new_idents)
seurat.shi.bc$cluster_level_1 <- Idents(seurat.shi.bc)

# Reporting
seurat_resolution_test(seurat.shi.bc, 
                       'cluster_level_1', 
                       0.4,
                       'pcw',
                       our_markers,
                       'shi_subset')

# Pull out differentially expressed markers
dEx.markers <- FindAllMarkers(seurat.shi.bc, only.pos = TRUE, 
                              min.pct = 0.25, logfc.threshold = 0.25)
top10 <- dEx.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

### captures count of cells for the ident with the fewest
maxcells  <- min(table(Idents(seurat.shi.bc)))
                 
### nested object subsetting with downsampling
DoHeatmap(subset(ScaleData(seurat.shi.bc), downsample = maxcells), features = top10$gene)

## Sub-clustering  ---------------------------------------------------------------------
for (REGION in c('MGE', 'LGE', 'CGE', 'Progenitor')) {
  
  seurat.shi.bc.obj <- subset(seurat.shi.bc, subset = cluster_level_1 %in% REGION)
  seurat.shi.bc.obj <- FindVariableFeatures(seurat.shi.bc.obj, selection.method = "vst", nfeatures = 2000)
  seurat.shi.bc.obj <- RunFastMNN(object.list = SplitObject(seurat.shi.bc.obj, split.by = "pcw"))
  seurat.shi.bc.obj <- RunUMAP(seurat.shi.bc.obj, reduction = "mnn", dims = 1:10)
  seurat.shi.bc.obj <- FindNeighbors(seurat.shi.bc.obj, reduction = "mnn", dims = 1:10)
  seurat.shi.bc.obj <- FindClusters(seurat.shi.bc.obj, resolution = c(0.3, 0.4, 0.5, 0.6))
  
  seurat.shi.bc.obj$cluster_level_2 <- Idents(seurat.shi.bc.obj)
  
  seurat_resolution_test(seurat.shi.bc.obj, 
                         'cluster_level_2', 
                         c(0.3, 0.4, 0.5, 0.6),
                         'pcw',
                         our_markers,
                         paste0('shi_',  tolower(REGION)))
  
  assign(paste0('seurat_', tolower(REGION)), seurat.shi.bc.obj)


}


#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------