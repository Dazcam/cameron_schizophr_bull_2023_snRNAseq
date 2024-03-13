#--------------------------------------------------------------------------------------
#
#    Shi 2021 snRNAseq data - Seurat 
#
#--------------------------------------------------------------------------------------

##  Info  -----------------------------------------------------------------------------

# 1. Run Seurat pipeline for shi data (Seurat 4.3.0)
# 2. Sub-cluster Seurat object on MGE, LGE, CGE, and Progenitor cells
# 3. Adds cluster levels 1 and 2 to main Seurat object
# 4. Save main Seurat and sub-cluster objects 

##  Load Packages  --------------------------------------------------------------------
if (!require("Require")) install.packages("Require")
Require::Require(c("tidyverse", "readxl", "data.table", "ggdendro", "Seurat", 
                   "SeuratWrappers", "cowplot", "scCustomize", "rmarkdown", "SeuratDisk")) 

##  Set Variables  --------------------------------------------------------------------
DATA_DIR <- '~/Desktop/fetal_brain_snRNAseq_GE_270922/resources/'
OUT_DIR <- '~/Desktop/fetal_brain_snRNAseq_GE_270922/results/'
R_DIR <- paste0(OUT_DIR, '01R_objects/')
SCRIPT_DIR <- '~/Desktop/fetal_brain_snRNAseq_GE_270922/workflow/scripts/'


##  Load scripts and functions   -------------------------------------------------------
source(paste0(SCRIPT_DIR, 'snRNAseq_GE_prep_shi_data_for_Seurat.R'))

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
new_idents <- c('MGE', 'CGE', 'LGE', 'Progenitor', 'Progenitor', 
                'LGE', 'Early_InN', 'LGE', 'Progenitor', 'MGE', 
                'Progenitor', 'Microglia')
seurat.shi.bc <- Rename_Clusters(seurat.shi.bc, new_idents)
seurat.shi.bc$cluster_level_1 <- Idents(seurat.shi.bc)

## Subset Shi data - with batch correction sub-clusters  ------------------------------
# Not sub-clustering Microglia as only 242 cells
for (REGION in c('MGE', 'CGE', 'LGE', 'Progenitor', 'Early_InN')) {
  
  cat('\nRunning Seurat sub-cluster analysis for:', REGION, '...\n\n')
  
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
  
}

## Shi data - join subcluster annotations to main Seurat object  ----------------------
seurat_shi_meta <- seurat.shi.bc@meta.data %>%
  mutate(cluster_level_2 = rep(NA, dim(seurat.shi.bc)[2])) %>%
  rownames_to_column('Cells')

seurat_lvl2_meta <- seurat_shi_meta %>% 
  left_join(seurat_MGE_meta %>% select(Cells, cluster_level_2), by = c('Cells'), suffix = c('', '.1')) %>%
  left_join(seurat_LGE_meta %>% select(Cells, cluster_level_2), by = c('Cells'), suffix = c('', '.2')) %>%
  left_join(seurat_CGE_meta %>% select(Cells, cluster_level_2), by = c('Cells'), suffix = c('', '.3')) %>%
  left_join(seurat_Early_InN_meta %>% select(Cells, cluster_level_2), by = c('Cells'), suffix = c('', '.4')) %>%
  left_join(seurat_Progenitor_meta %>% select(Cells, cluster_level_2), by = c('Cells'), suffix = c('', '.5')) %>%
  mutate(cluster_level_2 = coalesce(cluster_level_2.1, cluster_level_2.2, cluster_level_2.3, 
                                    cluster_level_2.4, cluster_level_2.5)) %>%
  mutate(cluster_level_2 = replace_na(cluster_level_2, 'Microglia')) %>%
  select(-cluster_level_2.1, -cluster_level_2.2, -cluster_level_2.3, -cluster_level_2.4, 
         -cluster_level_2.5)
  
# Add subcluster annotations to original 
seurat.shi.bc$cluster_level_2 <- seurat_lvl2_meta$cluster_level_2

# Sanity check
# identical(colnames(seurat.shi.bc), seurat_shi_meta$Cells)

# Downsample R object
table(seurat.shi.bc$cluster_level_1) # Lowest cnt MG at 242 cells
table(seurat.shi.bc$cluster_level_2) # Dwnsmpl both levels to 242 cells
seurat.shi.bc_dwnSmpl_lvl1  <- subset(seurat.shi.bc, downsample = min(table(seurat.shi.bc$cluster_level_1)))

seurat.shi.bc_dwnSmpl_lvl2 <- seurat.shi.bc
Idents(seurat.shi.bc_dwnSmpl_lvl2) <- seurat.shi.bc_dwnSmpl_lvl2$cluster_level_2 
seurat.shi.bc_dwnSmpl_lvl2  <- subset(seurat.shi.bc_dwnSmpl_lvl2, downsample = 242)
table(seurat.shi.bc_dwnSmpl_lvl2$cluster_level_2) # 32 out of 36 cell types 242 cells 

# Save R objects
saveRDS(object = seurat.shi.bc, paste0(R_DIR, 'seurat_shi_bc.rds'))
saveRDS(object = seurat.shi.bc_dwnSmpl_lvl1, paste0(R_DIR, 'seurat_shi_bc_dwnSmpl_lvl1.rds'))
saveRDS(object = seurat.shi.bc_dwnSmpl_lvl2, paste0(R_DIR, 'seurat_shi_bc_dwnSmpl_lvl2.rds'))

for (REGION in c('MGE', 'CGE', 'LGE', 'Progenitor', 'Early_InN')) {
  
  seurat_object <- get(paste0('seurat_', REGION))
  saveRDS(seurat_object, paste0(R_DIR, 'seurat_shi_bc_', REGION, '.rds'))
  
}

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
