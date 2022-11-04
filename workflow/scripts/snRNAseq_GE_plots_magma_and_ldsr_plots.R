# -------------------------------------------------------------------------------------
#
#    snRNAseq MAGMA Celltyping and LDSR plots 
#
# -------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Code for figures: clusters lvl 1 and lvl 2
#  Issues: Region IDs encoded differently in MAGMA and LSDR analysis
#          Much tweaking needed between ggplot and cowplot for FDR < 5% line

##  Load packages  --------------------------------------------------------------------
cat('\nLoading packages ... \n\n')
library(tidyverse)
library(ggsignif)
library(cowplot)
library(reshape2) # For melt

##  Initialise variables  -------------------------------------------------------------
MAGMA_REGIONS <- c('cer', 'hip', 'pfc', 'tha', 'wge')
LDSR_REGIONS <- c("Cer", "FC", "GE", "Hipp", "Thal")
GWAS <- c('SCZ', 'HEIGHT')
RESULTS_DIR <- '~/Desktop/fetal_brain_snRNAseq_GE_270922/results/'
MAGMA_DATA_DIR <- paste0(RESULTS_DIR, 'magma/')
MAGMA_COND_DATA_DIR <- paste0(RESULTS_DIR, 'magma_conditional/')
LDSC_DATA_DIR <- paste0(RESULTS_DIR, 'LDSR_part_herit/baseline_v1.2/')
FIG_DIR <- paste0(RESULTS_DIR, 'figures/')
GENE_WINDOW <- c('10UP_10DOWN', '10UP_35DOWN')


# MAGMA - prepare df
cat('\nPreparing MAGMA data ... \n')
for (LEVEL in c('1', '2')) {

  for (DISORDER in GWAS) {
    
    for (WINDOW in GENE_WINDOW) {
    
    MAGMA_DF <- read.table(paste0(MAGMA_DATA_DIR, 'snRNAseq_GE_', DISORDER, '.lvl_', 
                                        LEVEL, '.magma.', WINDOW, '.gsa.out'), header = FALSE) %>%
      janitor::row_to_names(row_number = 1) %>% 
      mutate(VARIABLE = gsub('\\.', '-', VARIABLE)) %>%
      mutate(MAGMA = -log10(as.numeric(P))) %>%
      select(VARIABLE, MAGMA) %>%
      dplyr::rename(Category = VARIABLE) %>% # Match LDSRs cell-type column
      mutate(across('Category', str_replace, 'Early', 'Early_InN'))
             
    assign(paste0('magma_GE_', DISORDER, '_lvl_', LEVEL, '_', WINDOW, '_df'), MAGMA_DF, envir = .GlobalEnv) 
    
    }
 
  }
  
}

# LDSR - prepare df
cat('\nPreparing LDSR data ... \n')
for (DISORDER in GWAS) {
  
  # LDSR                                         
  LDSR_DF <- read_tsv(paste0(LDSC_DATA_DIR, 'snRNAseq_LDSR_', DISORDER, '_baseline.v1.2_summary.tsv')) %>%
    mutate(LDSR = if_else(`Coefficient_z-score` > 0, -log10(pnorm(`Coefficient_z-score`, lower.tail = FALSE)), 0)) %>%
    select(Category, LDSR)
    
  assign(paste0('ldsr_', DISORDER, '_df'), LDSR_DF, envir = .GlobalEnv) 
    
}
  
# Plot MAGMA and LDSR barplot
cat('\nCreate plots ... \n')
for (LEVEL in c('1', '2')) {
  
  for (DISORDER in GWAS) {
    
    if (LEVEL == '1') {
      
      LDSR_DF <- get(paste0('ldsr_', DISORDER, '_df')) %>%
        filter(Category %in% c('LGE', 'MGE', 'CGE', 'Early_InN', 'Progenitor', 'Microglia'))
      PLOT_DF <- left_join(get(paste0('magma_GE_', DISORDER, '_lvl_', LEVEL, '_df')), 
                           LDSR_DF, by = 'Category') %>% melt()
      BF_CORR <- 0.05/6
      
    } else {
      
      LDSR_DF <- get(paste0('ldsr_', DISORDER, '_df')) %>%
        filter(!Category %in% c('LGE', 'MGE', 'CGE', 'Early_InN', 'Progenitor', 'Microglia'))
      PLOT_DF <- left_join(get(paste0('magma_GE_', DISORDER, '_lvl_', LEVEL, '_df')), 
                           LDSR_DF, by = 'Category') %>% melt()
      BF_CORR <- 0.05/30
      
    }
        
    MAGMA_LDSR_PLOT <- ggplot(data = PLOT_DF, aes(x = value, y = factor(Category, rev(levels(factor(Category)))), 
                                                  fill = variable, group = rev(variable))) +
    geom_bar(stat = "identity", color = 'black', position = "dodge") +
      geom_vline(xintercept=-log10(BF_CORR), linetype = "dashed", color = "black") +
      geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
      theme_bw() +
      ggtitle(DISORDER) +
      theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", size = 1),
            plot.title = element_text(hjust = 0.5, face = 'bold'),
            axis.title.x = element_text(colour = "#000000", size = 14),
            axis.title.y = element_text(colour = "#000000", size = 14),
            axis.text.x  = element_text(colour = "#000000", size = 13, vjust = 0.5),
            axis.text.y  = element_text(colour = "#000000", size = 13),
            legend.position = "none") +
      xlab(expression(-log[10](P))) +
      ylab('Cell type') +
      xlim(0, 11.5) 
      
    assign(paste0(DISORDER, '_magma_ldsr_lvl_', LEVEL, '_plot'), MAGMA_LDSR_PLOT, envir = .GlobalEnv) 
    
  }

}


plot_grid(SCZ_magma_ldsr_lvl_1_plot, HEIGHT_magma_ldsr_lvl_1_plot)
plot_grid(SCZ_magma_ldsr_lvl_2_plot, HEIGHT_magma_ldsr_lvl_2_plot)
plot_grid(magma_cond_GE_MSN_plot, magma_cond_GE_InN_plot)
plot_grid(SCZ_magma_ldsr_lvl_2_plot, HEIGHT_magma_ldsr_lvl_2_plot)

plot_grid(SCZ_magma_ldsr_lvl_2_plot, 
          plot_grid(magma_cond_GE_MSN_plot, 
                    magma_cond_GE_InN_plot, ncol = 1))


legend <- get_legend(
  # create some space to the left of the legend
  ggplot(data = PLOT_DF, aes(x = value, y = Category, fill = variable, group = rev(variable))) +
    geom_bar(stat = "identity", color = 'black', position = "dodge") +
    theme(legend.key.size = unit(1.5, 'cm'),
          legend.text = element_text(size = 14),
          legend.title = element_blank()) 
)

# MAGMA conditional - prepare df
cat('\nPreparing MAGMA conditional data ... \n')
for (CELL_TYPE in c('InN', 'MSN')) {
  
  MAGMA_DF <- read.table(paste0(MAGMA_COND_DATA_DIR, 'magma_all_sig_and_skene_condition_skene_', 
                                CELL_TYPE, '.gsa.out'), header = TRUE) %>%
    mutate(MAGMA = -log10(as.numeric(P))) %>%
    select(VARIABLE, MAGMA) %>%
    dplyr::rename(Category = VARIABLE) %>%
    filter(!grepl('skene', Category))
  
  MAGMA_PLOT <- ggplot(data = MAGMA_DF, aes(x = MAGMA, y = factor(Category, rev(levels(factor(Category)))), 
                                            fill = '#F8766D')) +
    geom_bar(stat = "identity", color = 'black', position = "dodge") +
    geom_vline(xintercept=-log10(0.05/30), linetype = "dashed", color = "black") +
    geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
    theme_bw() +
    ggtitle(CELL_TYPE) +
    theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", size = 1),
          plot.title = element_text(hjust = 0.5, face = 'bold'),
          axis.title.x = element_text(colour = "#000000", size = 14),
          axis.title.y = element_text(colour = "#000000", size = 14),
          axis.text.x  = element_text(colour = "#000000", size = 13, vjust = 0.5),
          axis.text.y  = element_text(colour = "#000000", size = 13),
          legend.position = "none") +
    xlab(expression(-log[10](P))) +
    ylab('Cell type') +
    xlim(0, 11.5) 
  
  
  assign(paste0('magma_cond_GE_', CELL_TYPE, '_df'), MAGMA_DF, envir = .GlobalEnv) 
  assign(paste0('magma_cond_GE_', CELL_TYPE, '_plot'), MAGMA_PLOT, envir = .GlobalEnv)
  
}
# MAGMA downsampled
# MAGMA - prepare df
cat('\nPreparing MAGMA data ... \n')
for (LEVEL in c('1')) {
  
  for (DISORDER in GWAS) {
    
    for (GENE_WINDOW in (c('10UP_10DOWN', '10UP_35DOWN'))) {
    
      MAGMA_DF <- read.table(paste0(MAGMA_DATA_DIR, 'snRNAseq_GE_', DISORDER, 
                                        '.shi_bc_dwnSmpl.lvl_', LEVEL, '.magma.', 
                                        GENE_WINDOW, '.gsa.out'), 
                                 header = FALSE) %>%
      janitor::row_to_names(row_number = 1) %>% 
      mutate(VARIABLE = gsub('\\.', '-', VARIABLE)) %>%
      mutate(MAGMA = -log10(as.numeric(P))) %>%
      select(VARIABLE, MAGMA) %>%
      dplyr::rename(Category = VARIABLE) %>% # Match LDSRs cell-type column
      mutate(across('Category', str_replace, 'Early', 'Early_InN'))
    
      MAGMA_PLOT <- ggplot(data = MAGMA_DF, aes(x = MAGMA, y = factor(Category, rev(levels(factor(Category)))), 
                                                   fill = '#F8766D')) +
      geom_bar(stat = "identity", color = 'black', position = "dodge") +
      geom_vline(xintercept=-log10(0.05/30), linetype = "dashed", color = "black") +
      geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
      theme_bw() +
      ggtitle(paste0(DISORDER, "_", GENE_WINDOW)) +
      theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", size = 1),
            plot.title = element_text(hjust = 0.5, face = 'bold'),
            axis.title.x = element_text(colour = "#000000", size = 14),
            axis.title.y = element_text(colour = "#000000", size = 14),
            axis.text.x  = element_text(colour = "#000000", size = 13, vjust = 0.5),
            axis.text.y  = element_text(colour = "#000000", size = 13),
            legend.position = "none") +
      xlab(expression(-log[10](P))) +
      ylab('Cell type') +
      xlim(0, 11.5) 
      
      assign(paste0('magma_dwnSmpl_GE_', DISORDER, '_lvl_', LEVEL, '_', GENE_WINDOW, '_df'), MAGMA_DF, envir = .GlobalEnv)   
      assign(paste0('magma_dwnSmpl_GE_', DISORDER, '_lvl_', LEVEL, '_', GENE_WINDOW, '_plot'), MAGMA_PLOT, envir = .GlobalEnv) 
    

    
    }
    
  }
  
}


plot_grid(magma_GE_SCZ_lvl_1_10UP_10DOWN_plot, magma_GE_SCZ_lvl_1_10UP_35DOWN_plot,
          magma_GE_HEIGHT_lvl_1_10UP_10DOWN_plot, magma_GE_HEIGHT_lvl_1_10UP_35DOWN_plot)

# Is there a correlation between cell number and MAGMA / LDSR results?
magma_GE_SCZ_lvl_2_df$cell_counts <- as.vector(table(seurat.shi.bc$cluster_level_2))

cor(magma_GE_SCZ_lvl_2_df[, unlist(lapply(magma_GE_SCZ_lvl_2_df, is.numeric))])

test <- left_join(magma_GE_SCZ_lvl_2_df, ldsr_SCZ_df %>%
                    filter(!Category %in% c('LGE', 'MGE', 'CGE', 'Early_InN', 'Progenitor', 'Microglia')))

cor.test(test$MAGMA, test$cell_counts)
cor.test(test$LDSR, test$cell_counts)
cor.test(test$LDSR, test$MAGMA)
cor(test[, unlist(lapply(test, is.numeric))])
plot(test$LDSR, test$cell_counts)
plot(test$MAGMA, test$cell_counts)
plot(test$MAGMA, test$LDSR)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------


