#--------------------------------------------------------------------------------------
#
#    Prep Shi 2021 snRNAseq data - create ctd and MAGMA/LDSR gene lists
#
#--------------------------------------------------------------------------------------

##  Info  -----------------------------------------------------------------------------

  # 1. Create ctd for level 1 and level 2 clusters
  # 2. Prep gene lists for MAGMA and LDSR (latter bed files)
  # 3. Prep gene lists for MAGMA and LDSR lvl 2 conditional analyses
  # 4. Prep gene lists for top 1000 genes, rather than top 10% specific genes, analyses

##  Load Packages  --------------------------------------------------------------------
if (!require("Require")) install.packages("Require")
Require::Require(c("tidyverse", "readxl", "data.table", "BiocManager", "ggdendro",
                   "Seurat", "biomaRt")) 
# BiocManager::install(c("EWCE", "AnnotationDbi", "org.Hs.eg.db", "scuttle", zellkonverter"))

##  Set Variables  --------------------------------------------------------------------
DATA_DIR <- '~/Desktop/fetal_brain_snRNAseq_GE_270922/resources/'
OUT_DIR <- '~/Desktop/fetal_brain_snRNAseq_GE_270922/results/'
R_DIR  <- paste0(OUT_DIR, '01R_objects/')
SHI_DIR <- paste0(DATA_DIR, 'raw_data/shi_et_al_2021/')
H5AD_DIR <- paste0(OUT_DIR, 'h5ad_objects/')
GENELIST_DIR <- paste0(OUT_DIR, '02GENE_LISTS/')
MAGMA_DIR <- paste0(GENELIST_DIR, 'MAGMA/')
LDSR_DIR <- paste0(GENELIST_DIR, 'LDSR/')
PUBLIC_DATA <- paste0(DATA_DIR, 'public_data/skene_bryois_InN_entrez_gene_list.tsv')
dir.create(paste0(DATA_DIR, 'refs/'))
mhc_coords <- c(6, 28510120, 33480577)

# Rm MHC genes -----------------------------------------------------------------------
# hg19 reference file was downloaded from MAGMA website
gene_coords_raw <- read_tsv(paste0(DATA_DIR, 'refs/NCBI37.3.gene.loc.txt'),
                                 col_names = FALSE, col_types = 'cciicc')
gene_coords <- read_tsv(paste0(DATA_DIR, 'refs/NCBI37.3.gene.loc.txt'),
                               col_names = FALSE, col_types = 'cciicc') %>%
  dplyr::rename(entrez = "X1", chr = "X2", start = 'X3', 
                end = 'X4', strand = "X5", hgnc = 'X6') %>%
  filter(!(((start > mhc_coords[2] & start < mhc_coords[3]) | (end > mhc_coords[2] & end < mhc_coords[3])) &
             chr == mhc_coords[1])) %>%
#  write_tsv(paste0(DATA_DIR, 'refs/NCBI37.3.MHCremoved.gene.loc.txt'), col_names = FALSE) %>%
  mutate(chr = paste0("chr", chr)) %>%
  dplyr::select(chr, start, end, entrez, hgnc) 

cat('\nMHC genes removed: ', nrow(gene_coords_raw) - nrow(gene_coords))

##  Create ctd object  ----------------------------------------------------------------
# Requires raw count gene matrix - needs to be genes x cell and annotation data
for (ROBJ in c("", "_dwnSmpl_lvl1", "_dwnSmpl_lvl2")) {
  
  cat('\n\nCreating CTD object for:', paste0('seurat.shi.bc', ROBJ))
  
  SEURAT_OBJ <- readRDS(paste0(R_DIR, 'seurat_shi_bc', ROBJ, '.rds'))
  RAW_COUNTS <- SEURAT_OBJ@assays$RNA@counts
  RAW_COUNTS_NO_MHC <- RAW_COUNTS[!(rownames(RAW_COUNTS) %in% mhc_genes_uniq), ]
  cat('\nMHC genes removed:', dim(RAW_COUNTS)[1] - dim(RAW_COUNTS_NO_MHC)[1])
  
  # Create annotations 
  annotations <- as.data.frame(cbind(as.vector(rownames(SEURAT_OBJ@meta.data)),
                                     as.vector(SEURAT_OBJ$cluster_level_1), 
                                     as.vector(SEURAT_OBJ$cluster_level_2)))
  colnames(annotations) <- c('cell_id', 'level1class', 'level2class')
  rownames(annotations) <- NULL
  annotLevels <- list(level1class = annotations$level1class, 
                      level2class = annotations$level2class)
  
  # Normalize - this is optional, was not used in the original EWCE publication
  cat('\nRunning CPM ... ', '\n')
  COUNTs_CPM <- edgeR::cpm(RAW_COUNTS_NO_MHC)
  
  cat('\nDropping uninformative genes cpm norm ... ', '\n\n')
  DROP_GENES_CPM <- EWCE::drop_uninformative_genes(
    exp = COUNTs_CPM, 
    input_species = "human",
    output_species = "human",
    level2annot = annotLevels$level2class) 
  
  cat('\nGene counts:',
      '\n\nRAW:', dim(RAW_COUNTS)[1],
      '\nRAW_NO_MHC:', dim(RAW_COUNTS_NO_MHC)[1],
      '\nCPM_DROP_GENES:', dim(DROP_GENES_CPM)[1])
  
  # Create object - saves ctd obj to folder
  cat('\nCreating CTD object ... \n\n')
  dir.create(R_DIR, showWarnings = FALSE)
  ctd <- EWCE::generate_celltype_data(exp = DROP_GENES_CPM, 
                                      annotLevels = annotLevels, 
                                      groupName = paste0('shi', ROBJ),
                                      savePath = R_DIR,
                                      numberOfBins = 10)

  
}

##  Write MAGMA/LDSR input files ------------------------------------------------------
for (CTD_EXT in c("", "_dwnSmpl_lvl1", "_dwnSmpl_lvl2")) {
  
  cat(paste0('\n\nRunning shi_bc', CTD_EXT, ' df ...\n\n'))
  
  # Levels only relevant to MAGMA as cluster names explicit bed file name for LDSR
  if (CTD_EXT == "_dwnSmpl_lvl1") {
    
    LEVELS <- 1 
    SUB_DIR <- 'shi_bc_dwnSmpl/'
    
  } else if (CTD_EXT == "_dwnSmpl_lvl2") {
    
    LEVELS <- 2
    SUB_DIR <- 'shi_bc_dwnSmpl/'
    
  } else {
    
    LEVELS <- c(1, 2)
    SUB_DIR <- 'shi_bc/'
    
  }
  
  for (LEVEL in LEVELS) { 
    
    cat('Running cluster level:', LEVEL, '...\n\n')
    load(paste0(R_DIR, 'ctd_shi', CTD_EXT, '.rda'))
    
    dir.create(paste0(GENELIST_DIR, SUB_DIR, 'MAGMA/'),  recursive = TRUE, showWarnings = FALSE)
    dir.create(paste0(GENELIST_DIR, SUB_DIR, 'LDSR/'), recursive = TRUE, showWarnings = FALSE)
  
    CELL_TYPES <- colnames(ctd[[LEVEL]]$specificity_quantiles)
    
    MAGMA <- as_tibble(as.matrix(ctd[[LEVEL]]$specificity_quantiles), rownames = 'hgnc') %>%
      inner_join(gene_coords) %>%
      pivot_longer(all_of(CELL_TYPES), names_to = 'cell_type', values_to = 'quantile') %>%
      filter(quantile == 10) %>%
      dplyr::select(cell_type, entrez) %>%
      with(., split(entrez, cell_type))
    
    for(i in names(MAGMA)) {
      
      cat(i, " ", paste(MAGMA[[i]], collapse = " "), "\n", 
          file = paste0(GENELIST_DIR, SUB_DIR, 'MAGMA/shi_top10_lvl_', LEVEL, '.txt')
          , sep = '', append = TRUE)
      
    }
    
    # Note that we can only set the windows for LDSR as input for magma is just a
    # gene list, not a genomic region. Windows for genes in gene list are derived in MAGMA.
    # You cannot set non-symmetrical gene windows in LDSR.
    # for (WINDOW in c('0UP_0DOWN', '10UP_10DOWN', '35UP_10DOWN', '100UP_100DOWN')) {
    #   
    #   if (WINDOW == '0UP_0DOWN') {
    #     
    #     UPSTREAM <- 0
    #     DOWNSTREAM <- 0
    #     
    #     } else if (WINDOW == '10UP_10DOWN') {
    #     
    #     UPSTREAM <- 10000
    #     DOWNSTREAM <- 10000
    #     
    #     } else if (WINDOW == '35UP_10DOWN') {
    #     
    #     UPSTREAM <- 35000
    #     DOWNSTREAM <- 10000
    #     
    #     } else {
    #     
    #     UPSTREAM <- 100000
    #     DOWNSTREAM <- 100000
    #     
    #     }
    #   
    #   LDSR <- as_tibble(as.matrix(ctd[[LEVEL]]$specificity_quantiles), rownames = 'hgnc') %>%
    #     inner_join(gene_coordinates) %>%
    #     pivot_longer(all_of(CELL_TYPES), names_to = 'cell_type', values_to = 'quantile') %>%
    #     filter(quantile == 10) %>%
    #     mutate(start = ifelse(start - UPSTREAM < 0, 0, start - UPSTREAM), end = end + DOWNSTREAM) %>%
    #     select(chr, start, end, entrez, cell_type) %>%
    #     group_by(cell_type) %>%
    #     group_walk(~ write_tsv(.x[,1:4], paste0(GENELIST_DIR, SUB_DIR, 'LDSR/', 
    #                                             .y$cell_type, '.', WINDOW, '.bed'), col_names = FALSE))
    #     
    # }
    # 
  }

}

# Prepare file for MAGMA conditional analyses - need to check which cells are sig. first
source(paste0('~/Desktop/fetal_brain_snRNAseq_GE_270922/workflow/scripts/snRNAseq_GE_prep_conditional_genesets.R'))


# Additional analysis to run top n genes rather than top 10%
#load(paste0(R_DIR, 'ctd_shi.rda'))
for (CTD_EXT in c("")) {
  
  for (LEVEL in c(1, 2)) { 
    
    UPSTREAM <- 100000
    DOWNSTREAM <- 100000
    WINDOW <- '100UP_100DOWN'
    SUB_DIR <- 'shi_bc/'
    dir.create(paste0(GENELIST_DIR, SUB_DIR, 'LDSR_top_1000_genes/'), recursive = TRUE, showWarnings = FALSE)
    
    CELL_TYPES <- colnames(ctd[[LEVEL]]$specificity_quantiles)
  
  MAGMA <- as_tibble(ctd[[LEVEL]]$specificity, rownames = 'hgnc') %>%
    inner_join(gene_coords) %>%
    tidyr::pivot_longer(all_of(CELL_TYPES), names_to = 'cell_type', values_to = 'specificity') %>%
    group_by(cell_type) %>%
    top_n(n = 1000, wt = specificity) %>%
    dplyr::select(cell_type, entrez) %>%
    with(., split(entrez, cell_type))
  
  for(i in names(MAGMA)) {
    
    cat(i, " ", paste(MAGMA[[i]], collapse = " "), "\n", 
        file = paste0(GENELIST_DIR, SUB_DIR, 'MAGMA/shi_top1000_lvl_', LEVEL, '.txt')
        , sep = '', append = TRUE)
    
  }
  
  LDSR <- as_tibble(as.matrix(ctd[[LEVEL]]$specificity), rownames = 'hgnc') %>%
    inner_join(gene_coords) %>%
    pivot_longer(all_of(CELL_TYPES), names_to = 'cell_type', values_to = 'specificity') %>%
    mutate(start = ifelse(start - UPSTREAM < 0, 0, start - UPSTREAM), end = end + DOWNSTREAM) %>%
    group_by(cell_type) %>%
    top_n(n = 1000, wt = specificity) %>%
    dplyr::select(chr, start, end, entrez, cell_type) %>%
    group_walk(~ write_tsv(.x[,1:4], paste0(GENELIST_DIR, SUB_DIR, 'LDSR_top_1000_genes/', 
                                            .y$cell_type, '.', WINDOW, '.bed'), col_names = FALSE))
  
  }
  
}

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

              