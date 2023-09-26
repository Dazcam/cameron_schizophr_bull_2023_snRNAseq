# -------------------------------------------------------------------------------------
#
#    snRNAseq create tables for paper
#
# -------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

# S1 = Specificity values for all level 1 cell types
# S2 = Specificity values for all level 2 cell types
# S3-6 = Coordinates for OCRs for each broad cell population
# S7 = all_cells_PGC3_SCZ_finemapped_SNP_peak_overlaps_ext250bp_with_peaks_and_anns 
# S8 = Distal SNP PCR – promoter coaccessibility (Nick's table)

## Set variables  --------------------------------------------------------------------
RESULTS_DIR <- '~/Desktop/fetal_brain_snRNAseq_GE_270922/results/'
CTD_DIR <- paste0(RESULTS_DIR, '01R_objects/')
PEAK_DIR <- "~/Desktop/fetal_brain_snATACseq_V3_010323/results/05PEAKS/"
FINEMAPPED_DIR <- paste0(PEAK_DIR, "finemapped_SNPs/")
LVL1_CELL_TYPES <- c("CGE", "LGE", "MGE", "Progenitor")

## Load packages  --------------------------------------------------------------------
library(tidyverse)

# Get specificity scores for tables S1 and S2  --------------------------------------- 
load(paste0(CTD_DIR, 'ctd_shi.rda'))

spec_lvl_1 <- as.data.frame(as.matrix(ctd[[1]]$specificity)) %>%
  as_tibble(rownames = 'Gene') %>%
  rename(IPC = Early_InN) %>%
  rename_with(~ gsub("GE", "GE-N", .x, fixed = TRUE))

spec_lvl_2 <- as.data.frame(as.matrix(ctd[[2]]$specificity)) %>%
  as_tibble(rownames = 'Gene') %>%
  rename_with(~ gsub("GE", "GE-N", .x, fixed = TRUE)) %>%
  rename_with(~ gsub("_", "-", .x, fixed = TRUE)) 

# Get coordinates of all OCRs for each broad cell population table S3 -----------------
for (CELL_TYPE in LVL1_CELL_TYPES) {

  peaks_df <- read_tsv(paste0(PEAK_DIR, CELL_TYPE,'.hg38.ext250bp.bed'), col_names = FALSE) %>%
    select(X1, X2, X3) %>%
    rename(chr = X1) %>%
    rename(start = X2) %>%
    rename(end = X3) 
  
  # %>%
  #   rename_with(~ paste0(CELL_TYPE, '_chr'), 'X1') %>%
  #   rename_with(~ paste0(CELL_TYPE, '_start'), 'X2') %>%
  #   rename_with(~ paste0(CELL_TYPE, '_end'), 'X3') 

  assign(paste0(tolower(CELL_TYPE), '_peaks'), peaks_df)
  
}

# Get SNP peak overlaps binary table with annotations ---------------------------------
snp_overlaps <- read_tsv(paste0(FINEMAPPED_DIR, 'all_cells_PGC3_SCZ_finemapped_SNP_peak_overlaps_ext250bp_with_peaks_and_anns.tsv')) %>%
  select(-index_snp, -starts_with('gene'), -width) %>%
  relocate(rsid, chr, hg38_base_position, finemap_posterior_probability, start, end)

# Generate and save excel file  -------------------------------------------------------
excel_list <- list(spec_lvl_1, spec_lvl_2, cge_peaks, lge_peaks, mge_peaks, progenitor_peaks, snp_overlaps)

# Add names for excel tabs
for (i in seq(1, length(excel_list))) {
  
  names(excel_list)[[i]] <- paste0('S', i)
  
}

# Save list items in seprate tabs 
openxlsx::write.xlsx(excel_list, "~/Desktop/Supplementary_tables_v11.xlsx")

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------