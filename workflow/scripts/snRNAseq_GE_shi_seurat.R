#--------------------------------------------------------------------------------------
#
#    Shi 2021 snRNAseq data - Seurat 
#
#--------------------------------------------------------------------------------------

##  Info  -----------------------------------------------------------------------------

# 1. Create Seurat object and QC markdown file for shi data - MGE cells only
# 2. Sub-clusters Seurat object on MGE, LGE, CGE, and Progenitor cells
# 3. Adds cluster levels 1 and 2 to main Seurat object
# 4. Save main Seurat objects and objects for sub-clusters
# 5. Converts Seurat objects to h5ad format for scDRS

##  Load Packages  --------------------------------------------------------------------
if (!require("Require")) install.packages("Require")
Require::Require(c("tidyverse", "readxl", "data.table", "ggdendro", "Seurat", 
                   "SeuratWrappers", "cowplot", "scCustomize", "rmarkdown", "SeuratDisk")) 
# scCustomize for contrast in cols
# BiocManager::install(c("EWCE", "AnnotationDbi", "org.Hs.eg.db", "scuttle", zellkonverter"))

##  Set Variables  --------------------------------------------------------------------
DATA_DIR <- '~/Desktop/fetal_brain_snRNAseq_GE_270922/resources/'
OUT_DIR <- '~/Desktop/fetal_brain_snRNAseq_GE_270922/results/'
R_DIR <- paste0(OUT_DIR, 'R_objects/')
H5AD_DIR <- paste0(OUT_DIR, 'h5ad_objects/')
SCRIPT_DIR <- '~/Desktop/fetal_brain_snRNAseq_GE_270922/workflow/scripts/'
MARKDOWN_FILE <- paste0(SCRIPT_DIR, 'snRNAseq_GE_seurat.Rmd')
REPORT_DIR <- paste0(OUT_DIR, 'rmarkdown_reports/')
REPORT_FILE <- 'snRNAseq_GE_seurat.html'

##  Load scripts and functions   -------------------------------------------------------
source(paste0(SCRIPT_DIR, 'snRNAseq_GE_prep_shi_data_for_Seurat.R'))
source(paste0(SCRIPT_DIR, 'snRNAseq_GE_functions.R'))
source(paste0(SCRIPT_DIR, 'snRNAseq_GE_marker_genes.R'))

## Subset Shi data - with batch correction  ---------------------------------------------
seurat.shi <- CreateSeuratObject(counts = shi_data, meta.data = shi_meta)
seurat.shi.bc <- subset(seurat.shi, subset = pcw %in% c("09", "12_01", "13", "16", "18_01", 
                                                           "12_02", "12_02.1")) 
seurat.shi.bc$pcw <- str_replace(seurat.shi.bc$pcw, "12_02.1", "12_02")
seurat.shi.bc$ClusterID <- str_replace(seurat.shi.bc$ClusterID, "progenitor", "Progenitor")
seurat.shi.bc <- subset(seurat.shi.bc, subset = ClusterID %in% c("MGE","LGE", "CGE", "OPC", "Microglia",
                                                                 "Endothelial", "Progenitor"))
seurat.shi.bc <- NormalizeData(seurat.shi.bc)
seurat.shi.bc <- FindVariableFeatures(seurat.shi.bc, selection.method = "vst", nfeatures = 2000)
seurat.shi.bc <- RunFastMNN(object.list = SplitObject(seurat.shi.bc, split.by = "pcw"))
seurat.shi.bc <- RunUMAP(seurat.shi.bc, reduction = "mnn", dims = 1:10)
seurat.shi.bc <- FindNeighbors(seurat.shi.bc, reduction = "mnn", dims = 1:10)
seurat.shi.bc <- FindClusters(seurat.shi.bc, resolution = 0.5)

# Rename clusters
# Note we lose OPC and Endothelial here
new_idents <- c('MGE', 'CGE', 'LGE', 'Progenitor', 'Progenitor', 
                'LGE', 'Early_InN', 'LGE', 'Progenitor', 'MGE', 
                'Progenitor', 'Microglia')
seurat.shi.bc <- Rename_Clusters(seurat.shi.bc, new_idents)
seurat.shi.bc$cluster_level_1 <- Idents(seurat.shi.bc)

# Reporting
seurat_resolution_test(seurat.shi.bc, 
                       'cluster_level_1', 
                       0.5,
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
shi_subset_heatmap <- DoHeatmap(subset(ScaleData(seurat.shi.bc), 
                                       downsample = maxcells), features = top10$gene) 


# Create markdown file
render(MARKDOWN_FILE, output_file = REPORT_FILE, output_dir = REPORT_DIR)


## Subset Shi data - with batch correction sub-clusters  ------------------------------
for (REGION in c('MGE', 'LGE', 'CGE', 'Progenitor')) {
  
  cat('\nRunning Seurat sub-cluster analysis for:', REGION, '...\n\n')
  
  if (REGION == 'MGE') {
    
    MARKERS <- mge_level2_markers
    
  } else if (REGION == 'LGE') {
    
    MARKERS <- lge_level2_markers
    
  } else if (REGION == 'CGE') {
    
    MARKERS <- cge_markers
    
  } else {
    
    MARKERS <- progenitor_level2_markers 
    
  }
  
  SEURAT_OBJ <- subset(seurat.shi.bc, subset = cluster_level_1 %in% REGION)
  SEURAT_OBJ <- NormalizeData(SEURAT_OBJ)
  SEURAT_OBJ  <- FindVariableFeatures(SEURAT_OBJ, selection.method = "vst", nfeatures = 2000)
  SEURAT_OBJ  <- RunFastMNN(object.list = SplitObject(SEURAT_OBJ, split.by = "pcw"))
  SEURAT_OBJ  <- RunUMAP(SEURAT_OBJ , reduction = "mnn", dims = 1:10)
  SEURAT_OBJ  <- FindNeighbors(SEURAT_OBJ , reduction = "mnn", dims = 1:10)
  SEURAT_OBJ  <- FindClusters(SEURAT_OBJ, resolution = 0.5)
  
  SEURAT_META <- SEURAT_OBJ@meta.data %>%
    unite(cluster_level_2, c("cluster_level_1", "RNA_snn_res.0.5"),
          remove = FALSE) %>%
    rownames_to_column('Cells')
  
  assign(paste0('seurat_', REGION), SEURAT_OBJ, .GlobalEnv)
  assign(paste0('seurat_', REGION, '_meta'), SEURAT_META, .GlobalEnv)
  
  # cat('\nGenerating plots ...\n')
  # seurat_resolution_test(SEURAT_OBJ, 
  #                        'ClusterID',
  #                        0.5,
  #                        'pcw',
  #                        MARKERS,
  #                        paste0('shi_', tolower(REGION)))
  
}

DimPlot_scCustom(seurat.shi.bc, reduction = 'umap')
# Need to generate MARKDOWN for the final sub-clusterings

## Shi data - join subcluster annotations to main Seurat object  ----------------------
seurat_shi_meta <- seurat.shi.bc@meta.data %>%
  mutate(cluster_level_2 = rep(NA, dim(seurat.shi.bc)[2])) %>%
  rownames_to_column('Cells')

seurat_lvl2_meta <- seurat_shi_meta %>% 
  left_join(seurat_MGE_meta %>% select(Cells, cluster_level_2), by = c('Cells'), suffix = c('', '.1')) %>%
  left_join(seurat_LGE_meta %>% select(Cells, cluster_level_2), by = c('Cells'), suffix = c('', '.2')) %>%
  left_join(seurat_CGE_meta %>% select(Cells, cluster_level_2), by = c('Cells'), suffix = c('', '.3')) %>%
  left_join(seurat_Progenitor_meta %>% select(Cells, cluster_level_2), by = c('Cells'), suffix = c('', '.4')) %>%
  mutate(cluster_level_2 = coalesce(cluster_level_2, cluster_level_2.1, cluster_level_2.2, cluster_level_2.3, 
                                  cluster_level_2.4)) %>%
  mutate(cluster_level_2 = replace_na(cluster_level_2, 'Other')) %>%
  select(-cluster_level_2.1, -cluster_level_2.2, -cluster_level_2.3, -cluster_level_2.4)
  
# Add subcluster annotations to original 
seurat.shi.bc$cluster_level_2 <- seurat_lvl2_meta$cluster_level_2

# Sanity check
# identical(colnames(seurat.shi.bc), seurat_shi_meta$Cells)

# Save R object
saveRDS(object = seurat.shi.bc, paste0(R_DIR, 'seurat_shi_bc.rds'))

# Convert to h5ad object for scDRS
SaveH5Seurat(seurat.shi.bc, filename = paste0(H5AD_DIR, 'shi_bc.h5Seurat'))
Convert(paste0(H5AD_DIR, 'shi_bc.h5Seurat'), dest = "h5ad")

# Save R object and h5ad clusters for subclusters
for (REGION in c('MGE', 'LGE', 'CGE', 'Progenitor')) {
  
  cat('\nGenerating h5ad file for:', REGION, '...\n\n')
  
  SEURAT_OBJ <- get(paste0('seurat_', REGION))
  
  saveRDS(SEURAT_OBJ, paste0(R_DIR, 'seurat_shi_bc_', REGION, '.rds'))
  
  # Convert to h5ad object for scDRS
  SaveH5Seurat(SEURAT_OBJ, filename = paste0(H5AD_DIR, 'shi_bc_', REGION, '.h5Seurat'))
  Convert(paste0(H5AD_DIR, 'shi_bc_', REGION, '.h5Seurat'), dest = "h5ad")
  
}

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------