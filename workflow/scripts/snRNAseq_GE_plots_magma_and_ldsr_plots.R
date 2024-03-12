# -------------------------------------------------------------------------------------
#
#    snRNAseq MAGMA and LDSR barcharts
#
# -------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Code for figures: clusters lvl 1 and lvl 2

##  Load packages  --------------------------------------------------------------------
cat('\nLoading packages ... \n\n')
library(tidyverse)
library(ggsignif)
library(cowplot)
library(reshape2) # For melt
library(Seurat) # For no legend

##  Initialise variables  -------------------------------------------------------------
GWAS <- c('SCZ', 'BPD', 'ASD', 'MDD', 'ADHD', 'HEIGHT')
RESULTS_DIR <- '~/Desktop/fetal_brain_snRNAseq_GE_270922/results/'
MAGMA_DATA_DIR <- paste0(RESULTS_DIR, '04MAGMA/')
MAGMA_COND_DATA_DIR <- paste0(MAGMA_DATA_DIR, 'magma_conditional/')
LDSC_DATA_DIR <- paste0(RESULTS_DIR, '05LDSR/part_herit/baseline_v1.2/')
FIG_DIR <- paste0(RESULTS_DIR, 'figures/')
GENE_WINDOW <- c('10UP_10DOWN', '35UP_10DOWN', '100UP_100DOWN')
COND_CELL_TYPES <- c('CGE_1', 'CGE_2', 'LGE_1', 'LGE_2', 'LGE_4', 'MGE_2', 'MGE_3')
LVL1_CELL_TYPES <- c("Microglia", "Progenitor", "IPC", "MGE-N", "LGE-N", 
                     "CGE-N")
LVL2_CELL_TYPES <- c("Progenitor-10", "Progenitor-9", "Progenitor-8", "Progenitor-7", 
                     "Progenitor-6", "Progenitor-5", "Progenitor-4", "Progenitor-3", 
                     "Progenitor-2", "Progenitor-1", "Progenitor-0", "Microglia", "IPC-5", "IPC-4", 
                     "IPC-3", "IPC-2", "IPC-1", "IPC-0", "MGE-N-5", "MGE-N-4", "MGE-N-3", 
                     "MGE-N-2", "MGE-N-1", "MGE-N-0", "LGE-N-7", "LGE-N-6", "LGE-N-5", 
                     "LGE-N-4", "LGE-N-3", "LGE-N-2", "LGE-N-1", "LGE-N-0", "CGE-N-3", 
                     "CGE-N-2", "CGE-N-1", "CGE-N-0")
LVL1_CORR <- 6
LVL2_CORR <- 36

my_theme <- function() {
  
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", size = 1),
          plot.title = element_text(hjust = 0.5, face = 'bold'),
          axis.title.x = element_text(colour = "#000000", size = 14),
          axis.title.y = element_text(colour = "#000000", size = 14),
          axis.text.x  = element_text(colour = "#000000", size = 13, vjust = 0.5),
          axis.text.y  = element_text(colour = "#000000", size = 13),
          legend.text = element_text(size = 13),
          legend.title = element_blank(),
          strip.text.y = element_blank()) 
  
}

## INDIVIDUAL FIGURES  ----------------------------------------------------------------

## MAIN L1 and L2 FIGS -----

# MAGMA - prepare df
cat('\nPreparing MAGMA data ... \n')
for (LEVEL in c('1', '2')) {

  for (DISORDER in GWAS) {
        
    CORR <- if (LEVEL == '1') 0.05 / LVL1_CORR else 0.05 / LVL2_CORR

    MAGMA_DF <- read.table(paste0(MAGMA_DATA_DIR, 'snRNAseq_GE_', DISORDER, '.shi_bc.lvl_', 
                                        LEVEL, '.magma.35UP_10DOWN.gsa.out'), header = FALSE) %>%
      janitor::row_to_names(row_number = 1) %>% 
      mutate(VARIABLE = gsub('\\.', '-', VARIABLE)) %>%
      mutate(MAGMA = -log10(as.numeric(P))) %>%
      dplyr::select(VARIABLE, MAGMA) %>%
      dplyr::rename(Category = VARIABLE) %>%    
      mutate(Category = str_replace(Category, 'GE', 'GE-N')) %>%
      mutate(Category = str_replace(Category, 'Early_InN', 'IPC')) %>%
      mutate(Category = str_replace(Category, '_', '-')) 
             
    MAGMA_SIG <- MAGMA_DF %>%
      filter(MAGMA > -log10(CORR))
    
    assign(paste0('magma_', DISORDER, '_lvl_', LEVEL, '_35UP_10DOWN_df'), MAGMA_DF, envir = .GlobalEnv) 
    assign(paste0('magma_', DISORDER, '_lvl_', LEVEL, '_35UP_10DOWN_sig_df'), MAGMA_SIG, envir = .GlobalEnv)
                    
 
  }
  
}
  
# LDSR - prepare df
cat('\nPreparing LDSR data ... \n')
for (DISORDER in GWAS) {
  
  for (LEVEL in c('1', '2')) {
    
    CORR <- if (LEVEL == '1') 0.05 / LVL1_CORR else 0.05 / LVL2_CORR
    
    for (WINDOW in GENE_WINDOW) {
      
      if (LEVEL == '1') {
        
        # LDSR                                         
        LDSR_FULL_DF <- read_tsv(paste0(LDSC_DATA_DIR, 'snRNAseq_LDSR_', DISORDER, '_baseline.v1.2_summary.tsv')) %>%
          mutate(LDSR = if_else(`Coefficient_z-score` > 0, -log10(pnorm(`Coefficient_z-score`, lower.tail = FALSE)), 0)) %>%
          separate(Category, into = c('Category', 'Window'), sep = '\\.') %>%
          filter(Category %in% c('LGE', 'MGE', 'CGE', 'Early_InN', 'Progenitor', 'Microglia')) %>%
          mutate(Category = str_replace(Category, 'GE', 'GE-N')) %>%
          mutate(Category = str_replace(Category, 'Early_InN', 'IPC')) %>%
          mutate(Category = str_replace(Category, '_', '-')) %>%
          filter(Window == (!!WINDOW)) 
        
      } else {
    
        LDSR_FULL_DF <- read_tsv(paste0(LDSC_DATA_DIR, 'snRNAseq_LDSR_', DISORDER, '_baseline.v1.2_summary.tsv')) %>%
          mutate(LDSR = if_else(`Coefficient_z-score` > 0, -log10(pnorm(`Coefficient_z-score`, lower.tail = FALSE)), 0)) %>%
          separate(Category, into=c('Category', 'Window'), sep = '\\.') %>%
          filter(Window == (!!WINDOW)) %>%
          mutate(Category = str_replace(Category, 'GE', 'GE-N')) %>%
          mutate(Category = str_replace(Category, 'Early_InN', 'IPC')) %>%
          mutate(Category = str_replace(Category, '_', '-')) %>%
          filter(!grepl("^CGE-N$|^LGE-N$|^MGE-N$|^IPC$|^Progenitor$", Category))
        
      }
        
      LDSR_DF <- LDSR_FULL_DF %>%
        dplyr::select(Category, LDSR)
      
      LDSR_SIG <- LDSR_DF %>%
        filter(LDSR > -log10(CORR))
      
      assign(paste0('ldsr_', DISORDER, '_lvl_', LEVEL, '_', WINDOW, '_df'), LDSR_DF, envir = .GlobalEnv) 
      assign(paste0('ldsr_', DISORDER, '_lvl_', LEVEL, '_', WINDOW, '_full_df'), LDSR_FULL_DF, envir = .GlobalEnv)
      assign(paste0('ldsr_', DISORDER, '_lvl_', LEVEL, '_', WINDOW, '_sig_df'), LDSR_SIG, envir = .GlobalEnv)
      
    }
    
  }
    
}

# Plot MAGMA and LDSR barplot 
cat('\nCreate plots ... \n')
for (LEVEL in c('1', '2')) {
  
  for (DISORDER in GWAS) {
    
    LEVELS <- if (LEVEL == '1') LVL1_CELL_TYPES else LVL2_CELL_TYPES
    CORR <- if (LEVEL == '1') 0.05 / LVL1_CORR else  0.05 / LVL2_CORR
    LIMIT <- if (LEVEL == '1') 6 else  8
    
    PLOT_DF <- left_join(get(paste0('magma_', DISORDER, '_lvl_', LEVEL, '_35UP_10DOWN_df')), 
                         get(paste0('ldsr_', DISORDER, '_lvl_', LEVEL, '_100UP_100DOWN_df')),
                         by = 'Category') %>%
      mutate(Category = factor(Category, levels = LEVELS)) %>%
      rowwise() %>% 
      mutate(MEAN = mean(c(MAGMA, LDSR))) %>%
      mutate(COL = case_when(
        MAGMA > -log10(CORR) & LDSR > -log10(CORR) ~ "Both",
        MAGMA > -log10(CORR) & LDSR <= -log10(CORR) ~ "MAGMA",
        MAGMA <= -log10(CORR) & LDSR > -log10(CORR) ~ "SLDSR",
        MAGMA < -log10(CORR) & LDSR < -log10(CORR) ~ "None")) %>%
      mutate(COL = factor(COL, levels = c('Both', 'MAGMA', 'SLDSR', 'None'))) %>%
      mutate(TYPE = str_extract(Category, "^.{2}")) %>%
      mutate(TYPE = factor(TYPE, levels = c('CG', 'LG', 'MG', 'IP', 'Mi', 'Pr')))
    
    MAGMA_LDSR_PLOT <- ggplot(data = PLOT_DF, aes(x = MEAN, y = Category, 
                                                  fill = COL)) +
    geom_bar(stat = "identity", color = 'black', position = "dodge", width = 0.8) +
      geom_vline(xintercept=-log10(CORR), linetype = "dashed", color = "black") +
      geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
      scale_fill_manual(values = c("Both" = "#00BA38", "MAGMA" =  "yellow", 
                                   "SLDSR" =  "#00B0F6", "None" =  "lightgrey"), 
                        drop = FALSE) +
      ggtitle(DISORDER) +
      xlim(0, LIMIT) +
      theme_bw() +
      my_theme() +
      xlab(expression(Mean -log[10](P))) +
      ylab('Cell type') 
    
    if (LEVEL == 2) {
      
      MAGMA_LDSR_PLOT <- MAGMA_LDSR_PLOT  +
        facet_grid(rows = vars(TYPE), scales = 'free', space = 'free')
      
    }
      
    assign(paste0(DISORDER, '_magma_ldsr_lvl_', LEVEL, '_plot'), MAGMA_LDSR_PLOT, envir = .GlobalEnv) 
    
  }

}

## Gene window comparisons  -----------------------------------------------------------
# Plot MAGMA and LDSR barplot - Omiotting for 
# cat('\nCreate plots ... \n')
# for (WINDOW in GENE_WINDOW) {
#   
#   DENOM <- 6
#   
#   PLOT_TITLE <- WINDOW %>%
#     str_replace('_', ' - ')
#   
#   PLOT_DF <- left_join(get(paste0('magma_SCZ_lvl_1_', WINDOW, '_df')), 
#                        get(paste0('ldsr_SCZ_lvl_1_', WINDOW, '_df')),
#                        by = 'Category') %>%
#     reshape2::melt()
#   
#   MAGMA_LDSR_PLOT <- ggplot(data = PLOT_DF, aes(x = value, y = factor(Category, rev(levels(factor(Category)))), 
#                                                 fill = variable, group = rev(variable))) +
#     geom_bar(stat = "identity", color = 'black', position = "dodge") +
#     geom_vline(xintercept=-log10(0.05/DENOM), linetype = "dashed", color = "black") +
#     geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
#     theme_bw() +
#     ggtitle(PLOT_TITLE) +
#     theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
#           panel.grid.major = element_blank(), 
#           panel.grid.minor = element_blank(),
#           panel.border = element_rect(colour = "black", size = 1),
#           plot.title = element_text(hjust = 0.5, face = 'bold'),
#           axis.title.x = element_text(colour = "#000000", size = 14),
#           axis.title.y = element_text(colour = "#000000", size = 14),
#           axis.text.x  = element_text(colour = "#000000", size = 13, vjust = 0.5),
#           axis.text.y  = element_text(colour = "#000000", size = 13),
#           legend.text = element_text(size = 13),
#           legend.title = element_blank()) +
#     xlab(expression(-log[10](P))) +
#     ylab('Cell type') +
#     xlim(0, 8) 
#   
#   assign(paste0('magma_ldsr_gene_window_', WINDOW, '_plot'), MAGMA_LDSR_PLOT, envir = .GlobalEnv) 
#   
# }

## Downsampled plots  -----------------------------------------------------------------
# MAGMA - prepare df
cat('\nPreparing MAGMA data ... \n')
for (LEVEL in c('1', '2')) {
  
  for (DISORDER in 'SCZ') {
    
    for (WINDOW in '35UP_10DOWN') {
    
      MAGMA_DF <- read.table(paste0(MAGMA_DATA_DIR, 'snRNAseq_GE_', DISORDER, 
                                        '.shi_bc_dwnSmpl.lvl_', LEVEL, '.magma.', 
                                        WINDOW, '.gsa.out'), 
                                 header = FALSE) %>%
        janitor::row_to_names(row_number = 1) %>% 
        mutate(VARIABLE = gsub('\\.', '-', VARIABLE)) %>%
        mutate(MAGMA = -log10(as.numeric(P))) %>%
        select(VARIABLE, MAGMA) %>%
        dplyr::rename(Category = VARIABLE) %>%
        filter(!str_detect(Category, "Other")) %>%    
        mutate(Category = str_replace(Category, 'GE', 'GE-N')) %>%
        mutate(Category = str_replace(Category, 'Early_InN', 'IPC')) %>%
        mutate(Category = str_replace(Category, '_', '-'))
      
      
      assign(paste0('magma_dwnSmpl_', DISORDER, '_lvl_', LEVEL, '_', WINDOW, '_df'), MAGMA_DF, envir = .GlobalEnv)   
    
    }
    
  }
  
}

# LDSR downsampled
# LDSR - prepare df
cat('\nPreparing LDSR data ... \n')
for (DISORDER in 'SCZ') {
  
  for (LEVEL in c('1', '2')) {
    
    for (WINDOW in '100UP_100DOWN') {
      
      if (LEVEL == '1') {
        
        # LDSR                                         
        LDSR_FULL_DF <- read_tsv(paste0(LDSC_DATA_DIR, 'shi_bc_dwnSmpl/snRNAseq_LDSR_', DISORDER, '_baseline.v1.2_summary.tsv')) %>%
          mutate(LDSR = if_else(`Coefficient_z-score` > 0, -log10(pnorm(`Coefficient_z-score`, lower.tail = FALSE)), 0)) %>%
          separate(Category, into = c('Category', 'Window'), sep = '\\.') %>%
          filter(Category %in% c('LGE', 'MGE', 'CGE', 'Early_InN', 'Progenitor', 'Microglia')) %>%
          mutate(Category = str_replace(Category, 'GE', 'GE-N')) %>%
          mutate(Category = str_replace(Category, 'Early_InN', 'IPC')) %>%
          mutate(Category = str_replace(Category, '_', '-')) %>%
          filter(Window == (!!WINDOW)) 
        
      } else {
        
        LDSR_FULL_DF <- read_tsv(paste0(LDSC_DATA_DIR, 'shi_bc_dwnSmpl/snRNAseq_LDSR_', DISORDER, '_baseline.v1.2_summary.tsv')) %>%
          mutate(LDSR = if_else(`Coefficient_z-score` > 0, -log10(pnorm(`Coefficient_z-score`, lower.tail = FALSE)), 0)) %>%
          separate(Category, into=c('Category', 'Window'), sep = '\\.') %>%
          mutate(Category = str_replace(Category, 'GE', 'GE-N')) %>%
          mutate(Category = str_replace(Category, 'Early_InN', 'IPC')) %>%
          mutate(Category = str_replace(Category, '_', '-')) %>%
          filter(Window == (!!WINDOW)) 
        
      }
      
      LDSR_DF <- LDSR_FULL_DF %>%
        select(Category, LDSR)
      
      assign(paste0('ldsr_dwnSmpl_', DISORDER, '_lvl_', LEVEL, '_', WINDOW, '_df'), LDSR_DF, envir = .GlobalEnv) 
      assign(paste0('ldsr_dwnSmpl_', DISORDER, '_lvl_', LEVEL, '_', WINDOW, '_full_df'), LDSR_FULL_DF, envir = .GlobalEnv)
      
    }
    
  }
  
}

# Plot
for (LEVEL in c('1', '2')) {
  

  LEVELS <- if (LEVEL == '1') LVL1_CELL_TYPES else LVL2_CELL_TYPES
  CORR <- if (LEVEL == '1') 0.05 / LVL1_CORR else 0.05 / LVL2_CORR
  
  cat('\nLEVEL:', LEVEL, '; CORR:', -log10(CORR), '\n')
  
  PLOT_DF <- left_join(get(paste0('magma_dwnSmpl_SCZ_lvl_', LEVEL, '_35UP_10DOWN_df')), 
                       get(paste0('ldsr_dwnSmpl_SCZ_lvl_', LEVEL, '_100UP_100DOWN_df')),
                       by = 'Category') %>%
    mutate(Category = factor(Category, levels = LEVELS)) %>%
    reshape2::melt()
  
  MAGMA_LDSR_PLOT <- ggplot(data = PLOT_DF, aes(x = value, y = Category, 
                                                fill = variable, group = rev(variable))) +
    geom_bar(stat = "identity", color = 'black', position = "dodge") +
    geom_vline(xintercept=-log10(CORR), linetype = "dashed", color = "black") +
    geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
    theme_bw() +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", size = 1),
          plot.title = element_text(hjust = 0.5, face = 'bold'),
          axis.title.x = element_text(colour = "#000000", size = 14),
          axis.title.y = element_text(colour = "#000000", size = 14),
          axis.text.x  = element_text(colour = "#000000", size = 13, vjust = 0.5),
          axis.text.y  = element_text(colour = "#000000", size = 13),
          legend.text = element_text(size = 13),
          legend.title = element_blank()) +
    xlab(expression(Mean -log[10](P))) +
    ylab('Cell type') +
    scale_x_continuous(breaks=seq(0, 10, 2))
  
  assign(paste0('magma_ldsr_dwnSmpl_SCZ_lvl_', LEVEL, '_plot'), MAGMA_LDSR_PLOT, envir = .GlobalEnv) 
  
}

## Top 1000 genes  ---------------------------------------------------------------------
# MAGMA - prepare df
cat('\nPreparing MAGMA data ... \n')
for (LEVEL in c('1', '2')) {
  
  for (DISORDER in 'SCZ') {
    
    for (WINDOW in c('35UP_10DOWN')) {
      
      MAGMA_DF <- read.table(paste0(MAGMA_DATA_DIR, 'snRNAseq_GE_', DISORDER, 
                                    '.shi_bc.top1000.lvl_', LEVEL, '.magma.', 
                                    WINDOW, '.gsa.out'), 
                             header = FALSE) %>%
        janitor::row_to_names(row_number = 1) %>% 
        mutate(VARIABLE = gsub('\\.', '-', VARIABLE)) %>%
        mutate(MAGMA = -log10(as.numeric(P))) %>%
        select(VARIABLE, MAGMA) %>%
        dplyr::rename(Category = VARIABLE) %>%
        mutate(Category = str_replace(Category, 'GE', 'GE-N')) %>%
        mutate(Category = str_replace(Category, 'Early_InN', 'IPC')) %>%
        mutate(Category = str_replace(Category, '_', '-'))
      
      assign(paste0('magma_top1000_', DISORDER, '_lvl_', LEVEL, '_', WINDOW, '_df'), MAGMA_DF, envir = .GlobalEnv)   

    }
    
  }
  
}


# LDSR - prepare df
cat('\nPreparing LDSR data ... \n')
for (DISORDER in 'SCZ') {
  
  for (LEVEL in c('1', '2')) {
    
    for (WINDOW in c('100UP_100DOWN')) {
      
      if (LEVEL == 1) {
        
        # LDSR                                         
        LDSR_FULL_DF <- read_tsv(paste0(LDSC_DATA_DIR, '/LDSR_top_1000_genes/snRNAseq_LDSR_', DISORDER, '_baseline.v1.2_summary.tsv')) %>%
          mutate(LDSR = if_else(`Coefficient_z-score` > 0, -log10(pnorm(`Coefficient_z-score`, lower.tail = FALSE)), 0)) %>%
          separate(Category, into = c('Category', 'Window'), sep = '\\.') %>%
          filter(Category %in% c('LGE', 'MGE', 'CGE', 'Early_InN', 'Progenitor', 'Microglia')) %>%
          dplyr::mutate(across(Category, str_replace, 'Early_InN', 'IPC')) %>%
          dplyr::mutate(across(Category, str_replace, 'GE', 'GE-N')) %>%
          filter(Window == (!!WINDOW)) 
        
      } else {
        
        LDSR_FULL_DF <- read_tsv(paste0(LDSC_DATA_DIR, '/LDSR_top_1000_genes/snRNAseq_LDSR_', DISORDER, '_baseline.v1.2_summary.tsv')) %>%
          mutate(LDSR = if_else(`Coefficient_z-score` > 0, -log10(pnorm(`Coefficient_z-score`, lower.tail = FALSE)), 0)) %>%
          separate(Category, into=c('Category', 'Window'), sep = '\\.') %>%
          mutate(Category = str_replace(Category, 'GE', 'GE-N')) %>%
          mutate(Category = str_replace(Category, 'Early_InN', 'IPC')) %>%
          filter(!Category %in% c('LGE', 'MGE', 'CGE', 'Early_InN', 'Progenitor')) %>%
          mutate(Category = str_replace(Category, '_', '-'))  %>%
          filter(Window == (!!WINDOW)) 
        
      }
      
      LDSR_DF <- LDSR_FULL_DF %>%
        select(Category, LDSR)
      
      
      assign(paste0('ldsr_top1000_', DISORDER, '_lvl_', LEVEL, '_', WINDOW, '_df'), LDSR_DF, envir = .GlobalEnv) 
      assign(paste0('ldsr_top1000_', DISORDER, '_lvl_', LEVEL, '_', WINDOW, '_full_df'), LDSR_FULL_DF, envir = .GlobalEnv)
      
    }
    
  }
  
}

# Plot
for (LEVEL in c('1', '2')) {
  
  LEVELS <- if (LEVEL == '1') LVL1_CELL_TYPES else LVL2_CELL_TYPES
  CORR <- if (LEVEL == '1') 0.05 / LVL1_CORR else 0.05 / LVL2_CORR
  
  cat('\nLEVEL:', LEVEL, '; CORR:', -log10(CORR), '\n')
  
  PLOT_DF <- left_join(get(paste0('magma_top1000_SCZ_lvl_', LEVEL, '_35UP_10DOWN_df')), 
                       get(paste0('ldsr_top1000_SCZ_lvl_', LEVEL, '_100UP_100DOWN_df')),
                       by = 'Category') %>%
    mutate(Category = factor(Category, levels = LEVELS)) %>%
    reshape2::melt()

  MAGMA_LDSR_PLOT <- ggplot(data = PLOT_DF, aes(x = value, y = Category, 
                                                fill = variable, group = rev(variable))) +
    geom_bar(stat = "identity", color = 'black', position = "dodge") +
    geom_vline(xintercept=-log10(CORR), linetype = "dashed", color = "black") +
    geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
    theme_bw() +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", size = 1),
          plot.title = element_text(hjust = 0.5, face = 'bold'),
          axis.title.x = element_text(colour = "#000000", size = 14),
          axis.title.y = element_text(colour = "#000000", size = 14),
          axis.text.x  = element_text(colour = "#000000", size = 13, vjust = 0.5),
          axis.text.y  = element_text(colour = "#000000", size = 13),
          legend.text = element_text(size = 13),
          legend.title = element_blank()) +
    xlab(expression(-log[10](P))) +
    ylab('Cell type') +
    xlim(0, 8) 
  
  assign(paste0('magma_ldsr_top1000_SCZ_lvl_', LEVEL, '_plot'), MAGMA_LDSR_PLOT, envir = .GlobalEnv) 
  
}


## Conditional analyses  --------------------------------------------------------------
# MAGMA conditional - prepare df
# cat('\nPreparing MAGMA conditional data ... \n')
# for (CELL_TYPE in COND_CELL_TYPES) {
#   
#   for (WINDOW in c('35UP_10DOWN')) { 
#     
#     MAGMA_DF <- read.table(paste0(MAGMA_COND_DATA_DIR, 'magma_all_sig_and_skene_condition_', 
#                                   CELL_TYPE, '.', WINDOW, '.gsa.out'), header = TRUE) %>%
#       mutate(MAGMA = -log10(as.numeric(P))) %>%
#       filter(!grepl(CELL_TYPE, VARIABLE)) %>%
#       filter(!grepl('skene', VARIABLE)) %>%
#       filter(!grepl('bryois', VARIABLE)) %>%
#       select(VARIABLE, MAGMA) %>%
#       dplyr::rename(Category = VARIABLE) %>%
#       dplyr::mutate(across(Category, str_replace, 'GE_', 'GE-N-')) 
#     
#     assign(paste0('magma_cond_GE_', CELL_TYPE, '_', WINDOW, '_df'), MAGMA_DF, envir = .GlobalEnv) 
#     
#   }
#   
# }
# 
# 
# # LDSR conditional - prepare df
# cat('\nPreparing LDSR conditional data ... \n')
# for (COND_DIR in c('LDSR_cond_int/', 'LDSR_cond_adult/')) {
#   
#   WINDOW <- '100UP_100DOWN'
#   
#   if (COND_DIR == 'LDSR_cond_int/') {
#     
#     for (CELL_TYPE in c('CGE_1', 'CGE_2', 'LGE_1', 'LGE_2', 'LGE_4', 'MGE_2', 'MGE_3')) {
#       
#       LDSR_DF <- read_tsv(paste0(LDSC_DATA_DIR, COND_DIR, 'snRNAseq_LDSR_SCZ_baseline.v1.2_summary.tsv')) %>%
#         mutate(LDSR = if_else(`Coefficient_z-score` > 0, -log10(pnorm(`Coefficient_z-score`, lower.tail = FALSE)), 0)) %>%
#         separate(Category, into = c('Category', 'Window'), sep = '\\.') %>%
#         separate(Category, into = c('Category', 'Condition'), sep = '_vs_') %>%
#         dplyr::mutate(across(Condition, str_replace, 'GE_', 'GE-N-')) %>%
#         filter(Category %in% CELL_TYPE) %>%
#         select(Condition, LDSR) %>%
#         dplyr::rename(Category = Condition)
#       
#       assign(paste0('ldsr_cond_GE_', CELL_TYPE, '_', WINDOW, '_df'), LDSR_DF, envir = .GlobalEnv) 
#       
#     }
#     
#   } else {
#     
#     for (CELL_TYPE in c('skene_InN', 'skene_MSN', 'bryois_GABA1', 'bryois_GABA2')) {
#       
#       LDSR_DF <- read_tsv(paste0(LDSC_DATA_DIR, COND_DIR, 'snRNAseq_LDSR_SCZ_baseline.v1.2_summary.tsv')) %>%
#         mutate(LDSR = if_else(`Coefficient_z-score` > 0, -log10(pnorm(`Coefficient_z-score`, lower.tail = FALSE)), 0)) %>%
#         separate(Category, into = c('Category', 'Window'), sep = '\\.') %>%
#         separate(Category, into = c('Category', 'Condition'), sep = '_vs_') %>%
#         dplyr::mutate(across(Category, str_replace, 'GE_', 'GE-N-')) %>%
#         filter(Condition %in% CELL_TYPE) %>%
#         dplyr::mutate(across(Condition, str_replace, 'skene_', 'Adult-mouse-')) %>%
#         dplyr::mutate(across(Condition, str_replace, 'bryois_', 'Adult-human-')) %>%
#         dplyr::select(Category, LDSR) 
#       
#       assign(paste0('ldsr_cond_GE_', CELL_TYPE, '_', WINDOW, '_df'), LDSR_DF, envir = .GlobalEnv) 
#       
#     }
#     
#     
#   }
#   
# }
# 
# # Plot
# for (CELL_TYPE in COND_CELL_TYPES) {
#   
#   if (grepl('skene', CELL_TYPE)) {
#     
#     CORR <- 0.05 / 7 
#     TITLE <- str_replace(CELL_TYPE, 'skene_', 'Adult mouse ')
#     
#   } else if (grepl('bryois', CELL_TYPE)) {
#     
#     CORR <- 0.05 / 7 
#     TITLE <- str_replace(CELL_TYPE, 'bryois_GABA', 'Adult human InN-')
#     
#   } else {
#     
#     CORR <- 0.05 / 6
#     TITLE <- str_replace(CELL_TYPE, 'GE_', 'GE-N-')
#     
#   }
#   
#   PLOT_DF <- left_join(get(paste0('magma_cond_GE_', CELL_TYPE, '_35UP_10DOWN_df')), 
#                        get(paste0('ldsr_cond_GE_', CELL_TYPE, '_100UP_100DOWN_df')),
#                        by = 'Category') %>%
#     reshape2::melt()
#   
#   MAGMA_LDSR_PLOT <- ggplot(data = PLOT_DF, aes(x = value, y = factor(Category, rev(levels(factor(Category)))), 
#                                                 fill = variable, group = rev(variable))) +
#     geom_bar(stat = "identity", color = 'black', position = "dodge") +
#     geom_vline(xintercept=-log10(CORR), linetype = "dashed", color = "black") +
#     geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
#     theme_bw() +
#     ggtitle(TITLE) +
#     theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
#           panel.grid.major = element_blank(), 
#           panel.grid.minor = element_blank(),
#           panel.border = element_rect(colour = "black", size = 1),
#           plot.title = element_text(hjust = 0.5, face = 'bold'),
#           axis.title.x = element_text(colour = "#000000", size = 14),
#           axis.title.y = element_text(colour = "#000000", size = 14),
#           axis.text.x  = element_text(colour = "#000000", size = 13, vjust = 0.5),
#           axis.text.y  = element_text(colour = "#000000", size = 13),
#           legend.text = element_text(size = 13),
#           legend.title = element_blank()) +
#     xlab(expression(-log[10](P))) +
#     ylab('Cell type') +
#     xlim(0, 8) 
#   
#   assign(paste0('magma_ldsr_cond_', CELL_TYPE, '_plot'), MAGMA_LDSR_PLOT, envir = .GlobalEnv) 
#   
# }

## Conditional analyses  NEW --------------------------------------------------------------
# MAGMA conditional - prepare df
cat('\nPreparing MAGMA conditional data ... \n')
for (CELL_TYPE in COND_CELL_TYPES) {
  
  names <- read.table(paste0(MAGMA_COND_DATA_DIR, 'magma_conditional_', 
                             CELL_TYPE, '.35UP_10DOWN.gsa.out'), header = TRUE) %>%
    slice(seq(2, nrow(.), 2)) %>%
    dplyr::select(VARIABLE) %>%
    slice_head(n = 17) %>%
    pull(VARIABLE)
  
  MAGMA_DF <- read.table(paste0(MAGMA_COND_DATA_DIR, 'magma_conditional_', 
                                CELL_TYPE, '.35UP_10DOWN.gsa.out'), header = TRUE) %>%
    mutate(MAGMA = -log10(as.numeric(P))) %>%
    slice_head(n = 34) %>%
    slice(seq(1, nrow(.), 2)) %>%
    mutate(VARIABLE = names) %>%
    dplyr::rename(Category = VARIABLE) %>%
    dplyr::mutate(Category = str_replace_all(Category, '_', '-')) %>%
    dplyr::select(Category, MAGMA) 
  
  assign(paste0('magma_cond_public_GE_', CELL_TYPE, '_df'), MAGMA_DF, envir = .GlobalEnv) 
  
}


# LDSR conditional - prepare df
cat('\nPreparing LDSR conditional data ... \n')
for (COND in c('int', 'public')) {
  
  if (COND == 'int') {
    
    for (CELL_TYPE in c('CGE_1', 'CGE_2', 'LGE_1', 'LGE_2', 'LGE_4', 'MGE_2', 'MGE_3')) {
      
      LDSR_DF <- read_tsv(paste0(LDSC_DATA_DIR, 'LDSR_cond_', COND, '/snRNAseq_LDSR_SCZ_baseline.v1.2_summary.tsv')) %>%
        mutate(LDSR = if_else(`Coefficient_z-score` > 0, -log10(pnorm(`Coefficient_z-score`, lower.tail = FALSE)), 0)) %>%
        separate(Category, into = c('Category', 'Window'), sep = '\\.') %>%
        separate(Category, into = c('Category', 'Condition'), sep = '_vs_') %>%
        dplyr::mutate(across(Condition, str_replace, 'GE_', 'GE-N-')) %>%
        filter(Category %in% CELL_TYPE) %>%
        dplyr::select(Condition, LDSR) %>%
        dplyr::rename(Category = Condition)
      
      assign(paste0('ldsr_cond_', COND, '_GE_', CELL_TYPE, '_df'), LDSR_DF, envir = .GlobalEnv) 
      
    }
    
  } else {
    
    for (CELL_TYPE in c('CGE_1', 'CGE_2', 'LGE_1', 'LGE_2', 'LGE_4', 'MGE_2', 'MGE_3')) {
      
      LDSR_DF <- read_tsv(paste0(LDSC_DATA_DIR, 'LDSR_cond_', COND, '/snRNAseq_LDSR_SCZ_baseline.v1.2_summary.tsv')) %>%
        mutate(LDSR = if_else(`Coefficient_z-score` > 0, -log10(pnorm(`Coefficient_z-score`, lower.tail = FALSE)), 0)) %>%
        separate(Category, into = c('Category', 'Window'), sep = '\\.') %>%
        separate(Category, into = c('Category', 'Condition'), sep = '_vs_') %>%
        filter(Category %in% CELL_TYPE) %>%
        dplyr::mutate(Condition = str_replace_all(Condition, '_', '-')) %>%
        dplyr::select(Condition, LDSR) %>%
        dplyr::rename(Category = Condition)
      
      assign(paste0('ldsr_cond_', COND, '_GE_', CELL_TYPE, '_df'), LDSR_DF, envir = .GlobalEnv) 
      
    }
    
    
  }
  
}

cond_plot_list <- list()

# Plot
for (CELL_TYPE in COND_CELL_TYPES) {
  
  for (COND in c('public')) {
    
    CORR <- if (COND == 'int') 0.05 / 6 else 0.05 / 17
    
    TITLE <- CELL_TYPE
    
    PLOT_DF <- left_join(get(paste0('magma_cond_', COND, '_GE_', CELL_TYPE, '_df')), 
                         get(paste0('ldsr_cond_', COND, '_GE_', CELL_TYPE, '_df')),
                         by = 'Category') %>%
      reshape2::melt()
    
    MAGMA_LDSR_PLOT <- ggplot(data = PLOT_DF, aes(x = value, y = factor(Category, rev(levels(factor(Category)))), 
                                                  fill = variable, group = rev(variable))) +
      geom_bar(stat = "identity", color = 'black', position = "dodge") +
      geom_vline(xintercept=-log10(CORR), linetype = "dashed", color = "black") +
      geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
      theme_bw() +
      ggtitle(TITLE) +
      theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", size = 1),
            plot.title = element_text(hjust = 0.5, face = 'bold'),
            axis.title.x = element_text(colour = "#000000", size = 14),
            axis.title.y = element_text(colour = "#000000", size = 14),
            axis.text.x  = element_text(colour = "#000000", size = 13, vjust = 0.5),
            axis.text.y  = element_text(colour = "#000000", size = 13),
            legend.text = element_text(size = 13),
            legend.title = element_blank()) +
      xlab(expression(-log[10](P))) +
      ylab('Cell type') +
      xlim(0, 12) 
    
    assign(paste0('magma_ldsr_cond_', COND, '_', CELL_TYPE, '_plot'), MAGMA_LDSR_PLOT, envir = .GlobalEnv) 
    cond_plot_list <- list(cond_plot_list, MAGMA_LDSR_PLOT)
    
  }
  
}

plot_grid(magma_ldsr_cond_public_CGE_1_plot, magma_ldsr_cond_public_CGE_2_plot,
          magma_ldsr_cond_public_LGE_1_plot, magma_ldsr_cond_public_LGE_2_plot,
          magma_ldsr_cond_public_LGE_4_plot, magma_ldsr_cond_public_MGE_2_plot,
          magma_ldsr_cond_public_MGE_3_plot, ncol = 3)


# Produce final plots for paper. ------------------------------------------------------
legend <- get_legend(SCZ_magma_ldsr_lvl_1_plot)

# Main plots - fig 2 - Lvl 1 - Gene set enrichment barplots
figure_2 <- plot_grid(SCZ_magma_ldsr_lvl_1_plot + NoLegend() + ggtitle(NULL), 
                      legend, ncol = 2, rel_widths = c(1, 0.3))

# Main plots - fig 3 - Lvl 2 - Gene set enrichment barplots
figure_3 <- plot_grid(SCZ_magma_ldsr_lvl_2_plot + NoLegend() + ggtitle(NULL), 
                      legend, ncol = 2, rel_widths = c(1, 0.3))

# Supp figs all GWAS
# Main plots - fig 2 - Lvl 1 - Gene set enrichment barplots
figure_S4_gwas <- plot_grid(SCZ_magma_ldsr_lvl_1_plot + NoLegend() + ggtitle('Schizophrenia'),
                            ADHD_magma_ldsr_lvl_1_plot + NoLegend() + ggtitle('ADHD'),
                            ASD_magma_ldsr_lvl_1_plot + NoLegend() + ggtitle('Autism'),
                            BPD_magma_ldsr_lvl_1_plot + NoLegend() + ggtitle('Bipolar Disorder'),
                            MDD_magma_ldsr_lvl_1_plot + NoLegend() + ggtitle('Major Depressive Disorder'), 
                            HEIGHT_magma_ldsr_lvl_1_plot + NoLegend() + ggtitle('Height'), 
                            ncol = 3)
figure_S4 <- ggdraw(figure_S4_gwas) +
  draw_plot(legend, .67, .45, .5, .4)

figure_S11_gwas <- plot_grid(SCZ_magma_ldsr_lvl_2_plot + NoLegend() + ggtitle('Schizophrenia'),
                             ADHD_magma_ldsr_lvl_2_plot + NoLegend() + ggtitle('ADHD'),
                             ASD_magma_ldsr_lvl_2_plot + NoLegend() + ggtitle('Autism'),
                             BPD_magma_ldsr_lvl_2_plot + NoLegend() + ggtitle('Bipolar Disorder'),
                             MDD_magma_ldsr_lvl_2_plot + NoLegend() + ggtitle('Major Depressive Disorder'), 
                             HEIGHT_magma_ldsr_lvl_2_plot + NoLegend() + ggtitle('Height'), 
                             ncol = 3)
figure_S11 <- ggdraw(figure_S11_gwas) +
  draw_plot(legend, .67, .45, .5, .3)

# Main plots - fig 2 - Lvl 2 - Gene set enrichment barplots
figure_4 <- plot_grid(SCZ_magma_ldsr_lvl_2_plot + NoLegend(), 
                      HEIGHT_magma_ldsr_lvl_2_plot+ NoLegend(), 
                      legend, ncol = 3, rel_widths = c(1, 1, 0.5), 
                      labels = c('A', 'B', ''), label_size = 20)

# Alternate gene window plots  
# gene_windows_plot <- plot_grid(magma_ldsr_gene_window_10UP_10DOWN_plot + NoLegend(), 
#                                magma_ldsr_gene_window_35UP_10DOWN_plot + NoLegend(),
#                                magma_ldsr_gene_window_100UP_100DOWN_plot + NoLegend(),
#                                legend, ncol = 4, rel_widths = c(1, 1, 1, 0.5),
#                                labels = c('A', 'B', 'C', ''), label_size = 20)

# Supplementary plots - S1 and S2 - dwnsampl and top1000 genes
lvl_1_supp_plot <- plot_grid(magma_ldsr_dwnSmpl_SCZ_lvl_1_plot + NoLegend(), 
                             magma_ldsr_top1000_SCZ_lvl_1_plot + NoLegend(), legend,
                             ncol = 3, labels = c('A', 'B', ''), label_size = 20)


lvl_2_supp_plot <- plot_grid(magma_ldsr_dwnSmpl_SCZ_lvl_2_plot + NoLegend(), 
                             magma_ldsr_top1000_SCZ_lvl_2_plot + NoLegend(), legend,
                             ncol = 3, labels = c('A', 'B', ''), label_size = 20)


# Supplementary plots - S7 - Internal conditional
cond_int_plot <- plot_grid(magma_ldsr_cond_CGE_1_plot + NoLegend(),
                           magma_ldsr_cond_CGE_2_plot + NoLegend(),
                           magma_ldsr_cond_LGE_1_plot + NoLegend(), 
                           magma_ldsr_cond_LGE_2_plot + NoLegend(),
                           magma_ldsr_cond_LGE_4_plot + NoLegend(),
                           magma_ldsr_cond_MGE_2_plot + NoLegend(),
                           magma_ldsr_cond_MGE_3_plot + NoLegend(),
                           ncol = 3, legend, labels = c('A', 'B', 'C', 'D', 'E',
                                                        'F', 'G', ''), 
                           label_size = 20)

# Supplementary plots - S8 - Adult conditional
cond_adult_plot <- plot_grid(magma_ldsr_cond_skene_InN_plot + NoLegend(), 
                             magma_ldsr_cond_skene_MSN_plot + NoLegend(),
                             magma_ldsr_cond_bryois_GABA1_plot + NoLegend(),
                             magma_ldsr_cond_bryois_GABA2_plot + NoLegend(), legend,
                             ncol = 3, labels = c('A', 'B', 'C', 'D', ''), label_size = 20)

# Produce tables?
# table_list <- c(ldsr_SCZ_lvl_1_100UP_100DOWN_full_df, ldsr_SCZ_lvl_2_100UP_100DOWN_full_df,
#                 ldsr_HEIGHT_lvl_1_100UP_100DOWN_full_df, ldsr_HEIGHT_lvl_2_100UP_100DOWN_full_df,
#                 )
#   
#   for (i in seq(1, length(excel_list))) {
#     
#     names(excel_list)[[i]] <- paste0('S', i)
#     
#   }
# openxlsx::write.xlsx(excel_list, "~/Desktop/Supplementary_tables_v11.xlsx")
# 




# Test code to split barchart by region
#  Prep data
cat('\nCreate plots ... \n')
for (LEVEL in c('1', '2')) {
  
  for (DISORDER in GWAS) {
    
    if (LEVEL == '1') {
      
      BF_CORR <- 0.05/6
      TEST_DF <- left_join(get(paste0('magma_', DISORDER, '_lvl_', LEVEL, '_35UP_10DOWN_df')), 
                           get(paste0('ldsr_', DISORDER, '_lvl_', LEVEL, '_100UP_100DOWN_df')),
                           by = 'Category') %>%
        reshape2::melt()
      
      TEST_PLOT <- ggplot(TEST_DF, aes(y = Category, x = value, fill = variable)) + 
        geom_bar(position = position_dodge2(reverse = TRUE), stat = "identity")
      
    } else {
      
      BF_CORR <- 0.05/29
      TEST_DF <- left_join(get(paste0('magma_', DISORDER, '_lvl_', LEVEL, '_35UP_10DOWN_df')), 
                           get(paste0('ldsr_', DISORDER, '_lvl_', LEVEL, '_100UP_100DOWN_df')),
                           by = 'Category') %>%
        reshape2::melt() %>% 
        mutate(Group = str_split(Category, "-", simplify = TRUE)[ , 1])
      TEST_DF$Group <- factor(TEST_DF$Group, levels = c("CGE", "LGE", "MGE", "Progenitor"))
      TEST_DF$Category <- factor(TEST_DF$Category,
                                 levels = c("CGE-N-3", "CGE-N-2", "CGE-N-1", "CGE-N-0", "LGE-N-7", 
                                            "LGE-N-6", "LGE-N-5", "LGE-N-4", "LGE-N-3", "LGE-N-2", 
                                            "LGE-N-1", "LGE-N-0", "MGE-N-5", "MGE-N-4", "MGE-N-3", 
                                            "MGE-N-2", "MGE-N-1", "MGE-N-0","Progenitor-10", 
                                            "Progenitor-9", "Progenitor-8", "Progenitor-7",
                                            "Progenitor-6", "Progenitor-5", "Progenitor-4", 
                                            "Progenitor-3", "Progenitor-2", "Progenitor-1", 
                                            "Progenitor-0"))
  
      TEST_PLOT <- ggplot(TEST_DF, aes(y = Category, x = value, fill = variable, group = Group)) + 
        geom_bar(position = position_dodge2(reverse = TRUE), stat = "identity")
      
    }
    
    # ggplot(data = PLOT_DF, aes(x = value, y = factor(Category, rev(levels(factor(Category)))), 
    #                            fill = variable, group = rev(variable))) +
    #   geom_bar(stat = "identity", color = 'black', position = "dodge") +
    #   geom_vline(xintercept=-log10(BF_CORR), linetype = "dashed", color = "black") +
    
    TEST_PLOT <- TEST_PLOT +
      theme_bw() +
      facet_grid(rows = vars(Group), scales = 'free_y', space = 'free') +
      geom_vline(xintercept=-log10(0.05/BF_CORR), linetype = "dashed", color = "black") +
      geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
      theme_bw() +
      ggtitle(paste0(DISORDER, '_', WINDOW)) +
      theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", linewidth = 1),
            plot.title = element_text(hjust = 0.5, face = 'bold'),
            axis.title.x = element_text(colour = "#000000", size = 14),
            axis.title.y = element_text(colour = "#000000", size = 14),
            axis.text.x  = element_text(colour = "#000000", size = 13, vjust = 0.5),
            axis.text.y  = element_text(colour = "#000000", size = 13),
            strip.background = element_blank(),
            strip.text.y = element_blank()) +
      xlab(expression(-log[10](P))) +
      ylab('Cell type') +
      xlim(0, 8) 
    
    assign(paste0(DISORDER, '_magma_ldsr_lvl_', LEVEL, '_TEST_plot'), TEST_PLOT, envir = .GlobalEnv) 
    
  }
  
}

# Level 1 - Gene set enrichment barplots
TEST_lvl_1 <- plot_grid(SCZ_magma_ldsr_lvl_1_TEST_plot + NoLegend(), 
                        HEIGHT_magma_ldsr_lvl_1_TEST_plot + NoLegend(), 
                        legend, ncol = 3, rel_widths = c(1, 1, 0.5),  
                        labels = c('A', 'B', ''), label_size = 20)

# Level 2 - Gene set enrichment barplots
TEST_lvl_2 <- plot_grid(SCZ_magma_ldsr_lvl_2_TEST_plot + NoLegend(), 
                        HEIGHT_magma_ldsr_lvl_2_TEST_plot + NoLegend(), 
                        legend, ncol = 3, rel_widths = c(1, 1, 0.5), 
                        labels = c('A', 'B', ''), label_size = 20)


## Diff exp genes  ---------------------------------------------------------------------
# MAGMA - prepare df
cat('\nPreparing MAGMA data ... \n')
magma_diffexp_tbl <- read.table(paste0(MAGMA_DATA_DIR, 
                              'magma_diffexp/snRNAseq_GE_magma_diffexp.35UP_10DOWN.gsa.out'), 
                       header = FALSE) %>%
  as_tibble() %>%
  janitor::row_to_names(row_number = 1) %>% 
  dplyr::rename(Category = VARIABLE) %>%
  mutate(Category = str_replace_all(Category, '_', '-')) %>%
  mutate(Category = str_replace(Category, 'early-InN', 'IPC')) %>%
  mutate(MAGMA = -log10(as.numeric(P))) %>%
  dplyr::select(Category, MAGMA)  
      
# LDSR - prepare df
ldsr_diffexp_tbl <- read_tsv(paste0(LDSC_DATA_DIR, 
                                    '/LDSR_DIFFEXP/snRNAseq_LDSR_SCZ_baseline.v1.2_summary.tsv')) %>%
  as_tibble() %>%
  mutate(LDSR = if_else(`Coefficient_z-score` > 0, -log10(pnorm(`Coefficient_z-score`, lower.tail = FALSE)), 0)) %>%
  separate(Category, into = c('Category', 'Window'), sep = '\\.') %>%
  mutate(Category = str_replace_all(Category, '_', '-')) %>%
  mutate(Category = str_replace(Category, 'early-InN', 'IPC')) %>%
  dplyr::select(Category, LDSR) 

# Plot
levels <- c("LGE-vs-CGE", "CGE-vs-IPC", "CGEandLGE-vs-IPC", "MGE-vs-IPC", "IPC-vs-Prog")
magma_diffexp_tbl %>%
  left_join(ldsr_diffexp_tbl) %>%
  mutate(Category = factor(Category, levels = levels))  %>% # order y-axis
  pivot_longer(c(MAGMA, LDSR)) %>%
  ggplot(aes(x = value, y = Category, 
                             fill = name, group = rev(name))) +
  geom_bar(stat = "identity", color = 'black', position = "dodge") +
  geom_vline(xintercept=-log10(0.05/4), linetype = "dashed", color = "black") +
  geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
  theme_bw() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1),
        plot.title = element_text(hjust = 0.5, face = 'bold'),
        axis.title.x = element_text(colour = "#000000", size = 14),
        axis.title.y = element_text(colour = "#000000", size = 14),
        axis.text.x  = element_text(colour = "#000000", size = 13, vjust = 0.5),
        axis.text.y  = element_text(colour = "#000000", size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_blank()) +
  xlab(expression(-log[10](P))) +
  ylab('Cell type') +
  xlim(0, 16) 
  

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------


