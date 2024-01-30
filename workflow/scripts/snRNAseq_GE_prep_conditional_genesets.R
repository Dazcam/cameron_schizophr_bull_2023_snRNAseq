#--------------------------------------------------------------------------------------
#
#    Prepare public gene sets for MAGMA conditional analyses
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Prepare publically available gene sets for conditional analyses
#  bryois_human_adult_top10pc.txt provided in private correspondance with author 
#  See: snRNAseq_GE_prep_enrich_test_files.R

## Set variables  ---------------------------------------------------------------------
require(tidyverse)
require(readxl)

## Set variables  ---------------------------------------------------------------------
options(timeout = max(100000, getOption("timeout"))) # Prevent timeout errors 
out_dir <- '~/Desktop/fetal_brain_snRNAseq_GE_270922/results/'
data_dir <- '~/Desktop/fetal_brain_snRNAseq_GE_270922/resources/public_data/'
genelist_dir <- paste0(out_dir, '02GENE_LISTS/')
cond_dir <- paste0(genelist_dir, 'shi_bc/MAGMA_CONDITIONAL/')
ldsr_cond_dir <- paste0(genelist_dir, 'shi_bc/LDSR/')
dir.create(ldsr_cond_dir,  recursive = TRUE, showWarnings = FALSE)
sig_cell_types <- c("CGE_1", "CGE_2","LGE_1", "LGE_2", "LGE_4", 
                    "MGE_2", "MGE_3")
upstream <- 100000
downstream <- 100000
window <- "100UP_100DOWN"

## Download data ----------------------------------------------------------------------
# cameron_url <- 'https://www.biologicalpsychiatryjournal.com/cms/10.1016/j.biopsych.2022.06.033/attachment/4083fa33-e040-446f-ab5c-95d0a5de71ec/mmc1.xlsx'
# skene_url <- 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-018-0129-5/MediaObjects/41588_2018_129_MOESM3_ESM.xlsx'

#download.file(cameron_url, paste0(data_dir, 'cameron2023_supp_tables.xlsx'), mode = "wb")
#download.file(skene_url, paste0(data_dir, 'skene2018_supp_table_S4.xlsx'), mode = "wb")

## Load data  --------------------------------------------------------------------------
load(paste0(R_DIR, 'ctd_shi.rda'))
fetal_tbl <- read_excel(paste0(data_dir, 'cameron2023_supp_tables.xlsx'), sheet = 'Table S7', skip = 2)
mouse_mat <- read_excel(paste0(data_dir, 'skene2018_supp_table_S4.xlsx')) %>%
  dplyr::rename(gene = `...1`) %>%
  dplyr::select('gene', 'interneurons',	'Medium Spiny Neuron',	'pyramidal CA1',	'pyramidal SS') %>%
  dplyr::rename('mouse_InN' ='interneurons',
                'mouse_MSN' = 'Medium Spiny Neuron', 
                'mouse_CA1' = 'pyramidal CA1',
                'mouse_SS' = 'pyramidal SS') %>%
  tibble::column_to_rownames(var = "gene") %>%
  as.matrix() 

mouse_quant <- apply(mouse_mat, 2, FUN = EWCE::bin_columns_into_quantiles, 
                     numberOfBins = 10) %>%
  as_tibble() %>%
  mutate(gene = rownames(mouse_mat)) %>%
  relocate(gene)
  
adult_tbl <- read_tsv(paste0(data_dir, 'bryois_human_adult_top10pc.txt'), col_names = FALSE) %>%
  t() %>% 
  as.data.frame() %>% 
  janitor::row_to_names(1) %>%
  as_tibble() %>%
  mutate_if(is.character, str_trim) %>%
  dplyr::rename('adult_exCA1' ='exCA1',
                'adult_exCA3' = 'exCA3', 
                'adult_exDG' = 'exDG',
                'adult_exPFC1' = 'exPFC1',
                'adult_exPFC2' = 'exPFC2',
                'adult_GABA1' = 'GABA1',
                'adult_GABA2' = 'GABA2')

magma_cond_list <- list()

for (cell_type in c('FC-ExN-2', 'FC-ExN-3', 'FC-ExN-4', 
                    'GE-InN-2', 'Hipp-ExN-3', 'Hipp-ExN-5')) {
  
  new_cell_type <- paste0('fetal_', str_replace_all(cell_type, '-', '_'))
  
  entrez_list <- fetal_tbl %>%
    dplyr::select(!!cell_type) %>%
    na.omit()  %>%
    inner_join(gene_coords, by = join_by(!!cell_type == hgnc))  %>%
    dplyr::select(!!cell_type, entrez) %>%
    dplyr::rename(!!new_cell_type := !!cell_type) %>%
    with(., split(entrez, new_cell_type))
  
  magma_cond_list <- c(magma_cond_list, entrez_list)
  
# Checks
# genes_tbl <- fetal_tbl %>%
#   dplyr::select(!!cell_type) %>%
#   na.omit() 
# 
# genes_tbl <- inner_join(genes_tbl, gene_coords, by = join_by(!!cell_type == hgnc)) 
# message(cell_type, ':before:', nrow(MAGMA), ':after:', nrow(MAGMA_2))
  
}

# Extract top 10% specificty genes from sig. mouse cell-types in Skene 2018 and mgi IDs to entrez
# Convert mouse gene IDs to human
# Use archived data to avoid error messages: https://support.bioconductor.org/p/9144001/
message('Converting gene IDs from mouse to human using BiomaRt ... ')
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

# Often thriows error due to server issues
# mart <- useEnsembl(biomart = "ensembl", mirror = "useast")
# human <- useDataset('hsapiens_gene_ensembl', mart)
# mouse <- useDataset('mmusculus_gene_ensembl', mart)

for (cell_type in c('mouse_InN', 'mouse_MSN', 'mouse_CA1', 'mouse_SS')) {
  
  # Pull out Q10 vector of genes for each sig, mouse cell type
  message(paste0('Extracting genes that have highest expression specificity in: ', cell_type))
  genes <- mouse_quant %>%
    dplyr::select(gene, !!cell_type) %>%
    dplyr::filter(.data[[cell_type]] == 10) %>%
    pull(gene) 
  
  message(paste0('Total genes before converting IDs: ', length(genes)))
  
  gene_lookup <- getLDS(attributes = c("mgi_symbol"), 
                       filters = "mgi_symbol", 
                       values = genes, 
                       mart = mouse, 
                       attributesL = c("hgnc_symbol", "entrezgene_id"), 
                       martL = human, 
                       uniqueRows = T,
                       ) %>% 
    as_tibble()  %>% 
    janitor::clean_names() %>%
    dplyr::rename(mgi = mgi_symbol,
                  !!cell_type := hgnc_symbol,
                  entrez = ncbi_gene_formerly_entrezgene_id) %>%
    unique() %>%
    drop_na() %>%
    mutate(entrez = as.character(entrez))
  
  entrez_genes <- gene_lookup %>%
    inner_join(gene_coords, by = join_by(!!cell_type == hgnc, entrez)) %>%
    unique() %>%
    drop_na() %>%
    dplyr::select(!!cell_type, entrez) %>%
    with(., split(entrez, cell_type))
  
  magma_cond_list <- c(magma_cond_list, entrez_genes)
  
  message(paste0('Total genes after converting IDs: ', length(entrez_genes[[cell_type]])))
  
}

for (cell_type in c('adult_exCA1', 'adult_exCA3', 'adult_exPFC1', 'adult_exDG', 
                    'adult_exPFC2', 'adult_GABA1', 'adult_GABA2')) {
  
  entrez_list <- adult_tbl %>%
    dplyr::select(!!cell_type) %>%
    na.omit()  %>%
    inner_join(gene_coords, by = join_by(!!cell_type == entrez)) %>%
    dplyr::select(!!cell_type, hgnc) %>%
    dplyr::rename(entrez = !!cell_type, 
                  !!cell_type := hgnc) %>%
    with(., split(entrez, cell_type))
  
  magma_cond_list <- c(magma_cond_list, entrez_list)
  
}

# Add Shi sig. cell types to MAGMA conditional -----------------------------------------
entrez_list <- as_tibble(as.matrix(ctd[[2]]$specificity_quantiles), rownames = 'hgnc') %>%
  inner_join(gene_coords) %>%
  dplyr::select(chr, start, end, entrez, hgnc, any_of(sig_cell_types)) %>%
  pivot_longer(all_of(sig_cell_types), names_to = 'cell_type', values_to = 'quantile') %>%
  filter(quantile == 10) %>%
  dplyr::select(cell_type, entrez) %>%
  with(., split(entrez, cell_type))

magma_cond_list <- c(magma_cond_list, entrez_list)

# Write gene lists to file
for(i in names(magma_cond_list)) {
  
  cat(i, " ", paste(magma_cond_list[[i]], collapse = " "), "\n", 
      file = paste0(cond_dir, 'snRNAseq_GE_conditional_gene_sets.txt')
      , sep = '', append = TRUE)
  
}


# Prepare files for LDSR conditional analyses -----------------------------------------
for (LINE in seq(1, 17, 1)) {
  
  cond_genelists <- readLines(paste0(cond_dir, 'snRNAseq_GE_conditional_gene_sets.txt'))
  genelist_name <- unlist(strsplit(cond_genelists[as.integer(LINE)], " "))[1]
  
  cat('Obtaining gene coords for', genelist_name, '...\n\n')
  
  df <- unlist(strsplit(cond_genelists[as.integer(LINE)], " ")) %>% 
    as_tibble() %>%
    janitor::row_to_names(row_number = 1) %>%
    dplyr::rename(entrez = 1) %>%
    inner_join(gene_coords) %>%
    dplyr::mutate(start = ifelse(start - upstream < 0, 0, start - upstream), end = end + downstream) %>%
    dplyr::select(chr, start, end, entrez) %>%
    write_tsv(paste0(ldsr_cond_dir, genelist_name, '.', window, '.bed'), col_names = FALSE)
      
}
  
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
