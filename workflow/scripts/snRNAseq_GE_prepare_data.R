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
SHI_DIR <- paste0(DATA_DIR, 'raw_data/shi_et_al_2021/')
OUT_DIR <- '~/Desktop/fetal_brain_snRNAseq_GE_270922/results/'
CTD_DIR <- paste0(OUT_DIR, 'ctd_objects/')
H5AD_DIR <- paste0(OUT_DIR, 'h5ad_objects/')
GENELIST_DIR <- paste0(OUT_DIR, 'gene_lists/')
MAGMA_DIR <- paste0(GENELIST_DIR, 'MAGMA/')
LDSR_DIR <- paste0(GENELIST_DIR, 'LDSR/')

##  Load Data  ------------------------------------------------------------------------
shi_data <- fread(paste0(SHI_DIR, "GSE135827_GE_mat_raw_count_with_week_info.txt"))
shi_meta <- read_excel(paste0(SHI_DIR, "science.abj6641_tables_s2_to_s9/science.abj6641_table_s2.xlsx"), 
                       col_names = TRUE, 
                       skip = 1) # Note added skip here to get rid of nonsense 1st line in excel sheet


##  Check if the cell orders are identical in shi_meta and shi_data  ------------------
shi_data[1:10, 1:10] # Note that fread added V1 as a col name for the genes column
shi_meta[1:10, 1:5]

# Pull out shi_data header so we don't have to deal with entire matrix
shi_data_cell_IDs_df <- t(head(shi_data, 1)) %>% 
  as.data.frame() %>%
  rownames_to_column("ID") %>%
  dplyr::slice(-1) %>% # Get rid of V1 row
  separate(ID, c("cell_ID", "pcw"), ".GW") %>%
  dplyr::select(-V1) %>%
  as_tibble() # get rid of V1 column

# Get the metadata cell IDs and remove trailing number 
shi_meta_cell_IDs_df <- shi_meta %>%
  separate(Cells, c("cell_ID", "cell_ID_number"), "-", remove = FALSE) 

# Note that the cell IDs in the metadata table do not match that give in the 
# colnames of the data matrix. The former have 1-11 assigned non uniformly
unique(shi_meta_cell_IDs_df$cell_ID_number)
table(shi_meta_cell_IDs_df$cell_ID_number) # Note these numbers are not uniform

# This number might be important - check if there are any duplicated cell IDs 
# in meta when that number is removed
sum(duplicated(shi_meta_cell_IDs_df$cell_ID)) 

# It's likely the numbers refer to different sequencing runs. Cells from 
# different sequencing runs can be given the same ID as the barcodes
# in each GEM are reused. We need to make sure the cell ID for each cell is unique
# Now check cell IDs in shi_data
sum(duplicated(shi_data_cell_IDs_df$cell_ID)) 
  
# Good sign that duplicate cells match - now check if the cell IDs in shi_meta and 
# shi_data are in the same order
identical(shi_meta_cell_IDs_df$cell_ID, shi_data_cell_IDs_df$cell_ID)

# Now check the format that CreateSeuratObject() needs the data in
#?CreateSeuratObject()

# Get metadata into correct format - retain clusterIDs and pcw info
shi_meta <- cbind(shi_meta, shi_data_cell_IDs_df)  %>%
  column_to_rownames(var = "Cells") %>%
  dplyr::select(`Major types`, pcw) %>%
  dplyr::rename(ClusterID = `Major types`)

# Now the count matrix
shi_data <- shi_data %>% 
  column_to_rownames(var = "V1")

# We've already extracted the pcw info and we know the cell orderings are identical
# so we can change the cell names in shi_data to exactly match that in the meta data
colnames(shi_data) <- shi_meta_cell_IDs_df$Cells

# Extract list of cell IDs to remove from matrix
cells_to_extract <- shi_meta %>%
  rownames_to_column(var = 'cell_type') %>%
  mutate(test = !(ClusterID  %in% c('Excitatory IPC', 
                                    'Thalamic neurons', 
                                    'Excitatory neuron'))) %>%
  pull(test)

# Remove columns annotated to 3 cell types we want to remove
shi_data_filt <- shi_data[, cells_to_extract]

# Extract unwanted cell types from metadata
shi_meta_filt <- shi_meta %>%
  rownames_to_column(var = 'cell_type') %>%
  filter(!grepl('Excitatory IPC|Thalamic neurons|Excitatory neuron', ClusterID)) 

# Might need to check this removal from two separate dfs is OK. And that cols
# removed from the mat correspond to rows removed from df.
# Using Seurat might be better. 


##  Create h5ad object  ---------------------------------------------------------------
dir.create(H5AD_DIR)
library(SingleCellExperiment)

# Prep data for sce generation
shi_cell_IDs <- shi_meta_filt$cell_type
shi_meta_for_sce <- shi_meta_filt[, 2:3]
rownames(shi_meta_for_sce) <- shi_cell_IDs

# Create sce object as input for zellkonverter
sce <- SingleCellExperiment(list(counts = shi_data_filt), colData = shi_meta_for_sce)
sce <- scuttle::addPerCellQC(sce)
sce <- addPerFeatureQC(sce)

# Create h5ad file for scDRS
zellkonverter::writeH5AD(sce, paste0(H5AD_DIR, 'shi2021_filt.h5ad'))


#seurat.shi <- CreateSeuratObject(counts = shi_data_filt, meta.data = shi_meta_for_seurat)
# 
# # Now that the data is loaded get rid of all the objects we don't need
# rm(shi_data, shi_data_cell_IDs_df, shi_meta, shi_meta_cell_IDs_df)
# 
# # Save Seurat object
# # Use this to derive the ctd objects in the gene specificity script
# saveRDS(seurat.shi, paste0(OUT_DIR, "seurat_shi.rds"))

##  Create ctd object  ----------------------------------------------------------------
# Requires raw count gene matrix - needs to be genes x cell and annotation data
# Create annotations 
annotations <- as.data.frame(cbind(rownames(shi_meta_filt),
                                   shi_meta_filt$ClusterID, 
                                   shi_meta_filt$ClusterID))
colnames(annotations) <- c('cell_id', 'level1class', 'level2class')
rownames(annotations) <- NULL
annotLevels <- list(level1class = annotations$level2class, 
                    level2class = annotations$level2class)
##  Create object  ------------------------------------------------------------------
dir.create(CTD_DIR)
ctd <- EWCE::generate_celltype_data(exp = shi_data_filt, 
                              annotLevels = annotLevels, 
                              groupName = 'shi',
                              savePath = CTD_DIR)

load(paste0(CTD_DIR, 'ctd_shi.rda'))

# Clean up
rm(list = ls(pattern = '^shi_*'))

##  Create specificity scores  ------------------------------------------------------
# Pull out raw mean expression scores from ctd
exp_lvl <- as.data.frame(as.matrix(ctd[[1]]$mean_exp)) %>% 
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


##  Write MAGMA/LDSR input files ------------------------------------------------------
# Filter out genes with expression below 1 (uninformative genes)
dir.create(MAGMA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LDSR_DIR, recursive = TRUE, showWarnings = FALSE)

MAGMA <- exp_specificity %>% 
  filter(Expr_sum_mean > 1) %>% 
  group_by(Lvl) %>% 
  top_n(., n_genes_to_keep, specificity) %>%
  select(Lvl, ENTREZ) %>%
  ungroup() %>%
  mutate(counts = rep(1:n_genes_to_keep, times = 7)) %>%
  pivot_wider(names_from = counts, values_from = ENTREZ) %>%
  write_tsv(paste0(MAGMA_DIR, "shi_top10.txt"), col_names = F)

LDSR <- exp_specificity %>% filter(Expr_sum_mean > 1) %>% 
  group_by(Lvl) %>% 
  top_n(., n_genes_to_keep, specificity) %>%
  select(chr, start, end, ENTREZ) %>%
  group_walk(~ write_tsv(.x[,1:4], paste0(LDSR_DIR, .y$Lvl, ".bed"), col_names = FALSE))


# Create gene coordinate file for LDSR
LDSR %>% 
  ungroup() %>%
  select(ENTREZ, chr, start, end) %>%
  write_tsv(paste0(LDSR_DIR, "LDSR_gene_coords.bed"))


#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

