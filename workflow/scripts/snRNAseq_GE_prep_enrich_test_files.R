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

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

