#--------------------------------------------------------------------------------------
#
#    Prepare the Shi 2021 snRNAseq data
#
#--------------------------------------------------------------------------------------

##  Load Packages  --------------------------------------------------------------------
library(tidyverse) 
library(EWCE) # To create ctd object
library(AnnotationDbi)
library(org.Hs.eg.db)
library(patchwork)
library(data.table)
library(ClusterProfiler)
library(readxl)


##  Set Variables  --------------------------------------------------------------------
DATA_DIR <- '~/Desktop/fetal_brain_snRNAseq_GE_270922/resources/raw_data/shi_et_al_2021/'
OUT_DIR <- '~/Desktop/fetal_brain_snRNAseq_GE_270922/results/ctd_objects/'
dir.create(OUT_DIR)

##  Load Data  ------------------------------------------------------------------------
shi_data <- fread(paste0(DATA_DIR, "GSE135827_GE_mat_raw_count_with_week_info.txt"))
shi_meta <- read_excel(paste0(DATA_DIR, "science.abj6641_tables_s2_to_s9/science.abj6641_table_s2.xlsx"), 
                       col_names = TRUE, 
                       skip = 1) # Note added skip here to get rid of nonsense 1st line in excel sheet


##  Check if the cell orders are identical in shi_meta and shi_data  ------------------
shi_data[1:10, 1:10] # Note that fread added V1 as a col name for the genes column
shi_meta[1:10, 1:5]

# Pull out header shi_data so we don't have to deal with entire matrix
shi_data_cell_IDs_df <- t(head(shi_data, 1)) %>% 
  as.data.frame() %>%
  rownames_to_column("ID") %>%
  slice(-1) %>% # Get rid of V1 row
  separate(ID, c("cell_ID", "pcw"), ".GW") %>%
  rename::select(-V1) # get rid of V1 column

# Get the metadata cell IDs and remove trailing number 
shi_meta_cell_IDs_df <- shi_meta %>%
  separate(Cells, c("cell_ID", "cell_ID_number"), "-", remove = FALSE) 

# Check if the number assigned to the cells in the shi_meta table are unique
unique(shi_meta_cell_IDs_df$cell_ID_number)

# This number might be important - check if there are any duplicated cell IDs 
# when that number is removed
sum(duplicated(shi_meta_cell_IDs_df$cell_ID)) 

# It's likely the numbers refer to different sequencing runs. Cells from 
# different sequencing runs can be given the same ID as the barcodes
# in each GEM are reused. We need to make sure the cell ID for each cell is unique
# Now check cell IDs in shi_data
sum(duplicated(shi_data_cell_IDs_df$cell_ID)) 
  
# Good sign that these match - now check if the cell IDs in shi_meta and shi_data are in the 
# same order
identical(shi_meta_cell_IDs_df$cell_ID, shi_data_cell_IDs_df$cell_ID)

# Now check the format that CreateSeuratObject() needs the data in
#?CreateSeuratObject()

# Get metadata into correct format
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

cells_to_extract <- shi_meta %>%
  rownames_to_column(var = 'cell_type') %>%
  mutate(test = !(ClusterID  %in% c('Excitatory IPC','Thalamic neurons','Excitatory neuron'))) %>%
  pull(test)

shi_meta_filt <- shi_meta %>%
  rownames_to_column(var = 'cell_type') %>%
  filter(!grepl('Excitatory IPC|Thalamic neurons|Excitatory neuron', ClusterID)) 

# Remove columns annotated to 3 cell types we want to remove
M <- shi_data[, cells_to_extract]

# # Create seurat object
# seurat.shi <- CreateSeuratObject(counts = shi_data, meta.data = shi_meta)
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
annotations <- as.data.frame(cbind(rownames(shi_meta), 
                             shi_meta$ClusterID, 
                             shi_meta$ClusterID))
colnames(annotations) <- c('cell_id', 'level1class', 'level2class')
rownames(annotations) <- NULL
annotLevels <- list(level1class = annotations$level2class, 
                    level2class = annotations$level2class)
##  Load ctd object  ------------------------------------------------------------------
ctd <- generate_celltype_data(exp = shi_data, 
                              annotLevels = annotLevels, 
                              groupName = 'shi',
                              savePath = OUT_DIR)

load("~/Desktop/fetal_brain_snRNAseq_GE_270922/results/CellTypeData_shi.rda")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
