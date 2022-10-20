#--------------------------------------------------------------------------------------
#
#    Prep Shi 2021 snRNAseq data - create ctd and MAGMA/LDSR gene lists
#
#--------------------------------------------------------------------------------------

##  Load Packages  --------------------------------------------------------------------

  # 1. Check shi data, remove 3 non-GE cell-types, munge for EWCE ctd creation
  # 2. Create ctd - may be adding fine grained cell types
  # 3. Prep gene lists for MAGMA and LDSR (inc. gene coord file)

##  Load Packages  --------------------------------------------------------------------
if (!require("Require")) install.packages("Require")
Require::Require(c("tidyverse", "readxl", "data.table", "BiocManager", "ggdendro",
                   "Seurat", "reticulate", "sceasy")) # Last 2 for AnnData objects
# BiocManager::install(c("EWCE", "AnnotationDbi", "org.Hs.eg.db", "scuttle", zellkonverter"))

##  Set Variables  --------------------------------------------------------------------
DATA_DIR <- '~/Desktop/fetal_brain_snRNAseq_GE_270922/resources/'
OUT_DIR <- '~/Desktop/fetal_brain_snRNAseq_GE_270922/results/'
R_DIR  <- paste0(OUT_DIR, 'R_objects/')
SHI_DIR <- paste0(DATA_DIR, 'raw_data/shi_et_al_2021/')
CTD_DIR <- paste0(OUT_DIR, 'ctd_objects/')
H5AD_DIR <- paste0(OUT_DIR, 'h5ad_objects/')
GENELIST_DIR <- paste0(OUT_DIR, 'gene_lists/')
MAGMA_DIR <- paste0(GENELIST_DIR, 'MAGMA/')
LDSR_DIR <- paste0(GENELIST_DIR, 'LDSR/')

##  Load Data  ------------------------------------------------------------------------
seurat.shi.bc <- readRDS(paste0(R_DIR, 'seurat_shi_bc.rds'))

##  Create ctd object  ----------------------------------------------------------------
# Requires raw count gene matrix - needs to be genes x cell and annotation data
# Create annotations 
annotations <- as.data.frame(cbind(as.vector(rownames(seurat.shi.bc@meta.data)),
                                   as.vector(seurat.shi.bc$cluster_level_1), 
                                   as.vector(seurat.shi.bc$cluster_level_2)))
colnames(annotations) <- c('cell_id', 'level1class', 'level2class')
rownames(annotations) <- NULL
annotLevels <- list(level1class = annotations$level1class, 
                    level2class = annotations$level2class)

##  Create object  ------------------------------------------------------------------
dir.create(CTD_DIR)
ctd <- EWCE::generate_celltype_data(exp = seurat.shi.bc@assays$RNA@counts, 
                              annotLevels = annotLevels, 
                              groupName = 'shi',
                              savePath = CTD_DIR)

load(paste0(CTD_DIR, 'ctd_shi.rda'))

##  Create specificity scores  ------------------------------------------------------
# Pull out raw mean expression scores from ctd

for (LEVEL in c(1, 2)) { 
  
  exp_lvl <- as.data.frame(as.matrix(ctd[[LEVEL]]$mean_exp)) %>% 
    rownames_to_column("Gene")
  
  # Gather data
  exp_lvl <- exp_lvl %>%
    gather(key = Lvl, value = Expr_sum_mean, -Gene) %>%
    as_tibble()
  
  # Load gene coordinates
  # Load hg19 gene coordinates and extend upstream and downstream coordinates by 100kb.
  # File downloaded from MAGMA website (https://ctg.cncr.nl/software/magma).
  # Filtered to remove extended MHC (chr6, 25Mb to 34Mb).
  gene_coordinates <- 
    read_tsv(paste0(DATA_DIR, 'refs/NCBI37.3.gene.loc.extendedMHCexcluded.txt'),
             col_names = FALSE, col_types = 'cciicc') %>%
    dplyr::rename(chr = "X2", ENTREZ = "X1", chr = "X2", 
                  start = 'X3', end = 'X4', HGNC = 'X6') %>%
  #  mutate(start = ifelse(X3 - 100000 < 0, 0, X3 - 100000), end = X4 + 100000) %>%
    dplyr::select(chr, start, end, ENTREZ) %>% 
    mutate(chr = paste0("chr",chr))
  
  # Write dictionary for cell type names
  dic_lvl <- dplyr::select(exp_lvl, Lvl) %>% 
    base::unique() %>% 
    mutate(makenames = make.names(Lvl))
  
  # Scale each cell type to the same total number of molecules (1M)
  exp_scaled <- exp_lvl %>% 
    group_by(Lvl) %>% 
    mutate(Expr_sum_mean = Expr_sum_mean * 1e6 / sum(Expr_sum_mean))
  
  # Specificity Calculation
  exp_specificity <- exp_scaled %>% group_by(Gene) %>% 
    mutate(specificity = Expr_sum_mean / sum(Expr_sum_mean))
  
  # Sanity check
  # exp_specificity %>%
  #   ungroup() %>%
  #   filter(Lvl == 'CGE') %>%
  #   summarise(across(where(is.numeric), sum))
  
  
  # Only keep MAGMA genes 
  entrez2symbol <- AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egSYMBOL2EG) %>% 
    dplyr::rename(Gene = "symbol", ENTREZ = "gene_id")
  exp_specificity <- inner_join(exp_specificity, entrez2symbol, by = "Gene") 
  exp_specificity <- inner_join(exp_specificity, gene_coordinates, by = "ENTREZ") 
  
  
  # Get number of genes that represent 10% of the dataset
  n_genes <- length(unique(exp_specificity$ENTREZ))
  n_genes_to_keep <- (n_genes * 0.1) %>% round()
  
  assign(paste0('exp_specificity_lvl_', LEVEL), exp_specificity, .GlobalEnv)
  assign(paste0('n_genes_to_keep_lvl_', LEVEL), n_genes_to_keep, .GlobalEnv)

}

ggplot(exp_specificity_lvl_1, aes(x = -log10(specificity), colour = Lvl)) + 
  geom_density()

##  Write MAGMA/LDSR input files ------------------------------------------------------
# Filter out genes with expression below 1 (uninformative genes)
dir.create(MAGMA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LDSR_DIR, recursive = TRUE, showWarnings = FALSE)

for (LEVEL in c(1, 2)) { 
  
  EXP_SPECIFICITY_DF <- get(paste0('exp_specificity_lvl_', LEVEL))
  N_GENES_TO_KEEP <- get(paste0('n_genes_to_keep_lvl_', 1))
  CELL_TYPE_CNT <- length(unique(EXP_SPECIFICITY_DF$Lvl))
  
  cat('\n\nRunning cluster level:', LEVEL, '...\n')
  cat(N_GENES_TO_KEEP, 'genes per cluster ...\n' )
  cat('Total cell type count:', CELL_TYPE_CNT, '...\n' )

  MAGMA <- EXP_SPECIFICITY_DF %>% 
    filter(Expr_sum_mean > 1) %>% 
    group_by(Lvl) %>% 
    top_n(., N_GENES_TO_KEEP, specificity) %>%
    select(Lvl, ENTREZ) %>%
    ungroup() %>%
    mutate(counts = rep(1:N_GENES_TO_KEEP, times = CELL_TYPE_CNT)) %>%
    pivot_wider(names_from = counts, values_from = ENTREZ) %>%
    write_tsv(paste0(MAGMA_DIR, 'shi_top10_lvl_', LEVEL, '.txt'), col_names = F)
  
  LDSR <- EXP_SPECIFICITY_DF %>% filter(Expr_sum_mean > 1) %>% 
    group_by(Lvl) %>% 
    top_n(., N_GENES_TO_KEEP, specificity) %>%
    select(chr, start, end, ENTREZ) %>%
    group_walk(~ write_tsv(.x[,1:4], paste0(LDSR_DIR, .y$Lvl, ".bed"), col_names = FALSE))

}

# Don't need this if using bed files for gene sets
# Create gene coordinate file for LDSR
# LDSR %>% 
#   ungroup() %>%
#   select(ENTREZ, chr, start, end) %>%
#   write_tsv(paste0(LDSR_DIR, "LDSR_gene_coords.bed"))


#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------


# + **GAD1** - InN
# + **GAD2** - InN
# + **SLC32A1** - InN
# + **GLI3** - RG
# + **SLC17A7** - ExN
# + **TNC** - RG
# + **PROX1** - CGE
# + **SCGN** - CGE
# + **LHX6** - MGE
# + **NXPH1** - MGE
# + **MEIS2** -LGE
# + **ZFHX3** - LGE
# + **SPI1** - MG
# + **LHX8** - MGE
# + **ISL1** - LGE

# This is code to remove cells from raw GeX matrix and metadat files provided by Shi
# If we are using Seurat this is probs not required. Filtering in Seurat is a much 
# Cleaner way to filter the cells. 
# Extract list of cell IDs to remove from matrix
# cells_to_extract <- shi_meta %>%
#   rownames_to_column(var = 'cell_type') %>%
#   mutate(test = !(ClusterID  %in% c('Excitatory IPC', 
#                                     'Thalamic neurons', 
#                                     'Excitatory neuron'))) %>%
#   pull(test)
# 
# # Remove columns annotated to 3 cell types we want to remove
# shi_data_filt <- shi_data[, cells_to_extract]
# 
# # Extract unwanted cell types from metadata
# shi_meta_filt <- shi_meta %>%
#   rownames_to_column(var = 'cell_type') %>%
#   filter(!grepl('Excitatory IPC|Thalamic neurons|Excitatory neuron', ClusterID)) 



##  Create h5ad object  ---------------------------------------------------------------
# Zellkonverter causes an error when running calculation in Scanpy. Opted to create
# Object seprately in python but undecided as yet what the best apporach is
# dir.create(H5AD_DIR)
# library(SingleCellExperiment)
# 
# # Prep data for sce generation
# shi_cell_IDs <- shi_meta_filt$cell_type
# shi_meta_for_sce <- shi_meta_filt[, 2:3]
# rownames(shi_meta_for_sce) <- shi_cell_IDs
# 
# # Create sce object as input for zellkonverter
# sce <- SingleCellExperiment(list(counts = shi_data_filt), colData = shi_meta_for_sce)
# sce <- scuttle::addPerCellQC(sce)
# sce <- addPerFeatureQC(sce)
# 
# # Create h5ad file for scDRS
# zellkonverter::writeH5AD(sce, paste0(H5AD_DIR, 'shi2021_filt.h5ad'))


# This is code for visualising scDRS results in Seurat and should be moved elsewhere
# DATA_DIR = '/Users/darren/Desktop/fetal_brain_snRNAseq_GE_270922/results/'
# scDRS_DIR = paste0(DATA_DIR, 'scDRS/')
# score <- read_tsv(paste0(scDRS_DIR, 'SCZ.full_score.gz'))
# seurat.shi$scz <- score$norm_score
# my_palette <- CustomPalette(low = "blue", high = "red", mid = 'white', k = 50)
# FeaturePlot(object = seurat.shi, features = "scz", cols = my_palette)

# Code for uploading shi UMAP mappings. 3d coordinates provided???
# seurat.shi@reductions[["umap"]]@cell.embeddings[,"UMAP_1"] <- shi_meta$`UMAP-X`
# seurat.shi@reductions[["umap"]]@cell.embeddings[,"UMAP_2"] <- shi_meta$`UMAP-Y`
# DimPlot(seurat.shi, reduction = "umap", group.by = 'pcw',
#         label = TRUE, label.size = 5,
#         pt.size = 0.1, repel = TRUE) + ggtitle(NULL) 

# # Now that the data is loaded get rid of all the objects we don't need
# rm(shi_data, shi_data_cell_IDs_df, shi_meta, shi_meta_cell_IDs_df)
# 
# # Save Seurat object
# # Use this to derive the ctd objects in the gene specificity script

# saveRDS(seurat.shi, paste0(OUT_DIR, "seurat_shi.rds"))
