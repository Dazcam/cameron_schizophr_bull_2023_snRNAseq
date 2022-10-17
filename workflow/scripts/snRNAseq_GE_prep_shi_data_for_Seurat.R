#--------------------------------------------------------------------------------------
#
#    Shi 2021 snRNAseq data - perpare shi data for seurat
#
#--------------------------------------------------------------------------------------

##  Info  --------------------------------------------------------------------

# 1. Prepare Shi data for seurat

##  Set Variables  --------------------------------------------------------------------
RES_DIR <- '~/Desktop/fetal_brain_snRNAseq_GE_270922/resources/' 
SHI_DIR <- paste0(RES_DIR, 'raw_data/shi_et_al_2021/')

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

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------