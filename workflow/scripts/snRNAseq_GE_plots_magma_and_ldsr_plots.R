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
library(Seurat) # For no legend

##  Initialise variables  -------------------------------------------------------------
MAGMA_REGIONS <- c('cer', 'hip', 'pfc', 'tha', 'wge')
LDSR_REGIONS <- c("Cer", "FC", "GE", "Hipp", "Thal")
GWAS <- c('SCZ', 'HEIGHT')
RESULTS_DIR <- '~/Desktop/fetal_brain_snRNAseq_GE_270922/results/'
MAGMA_DATA_DIR <- paste0(RESULTS_DIR, 'magma/')
MAGMA_COND_DATA_DIR <- paste0(RESULTS_DIR, 'magma_conditional/')
LDSC_DATA_DIR <- paste0(RESULTS_DIR, 'LDSR_part_herit/baseline_v1.2/')
FIG_DIR <- paste0(RESULTS_DIR, 'figures/')
GENE_WINDOW <- c('0UP_0DOWN', '10UP_10DOWN', '35UP_10DOWN', '100UP_100DOWN')
COND_CELL_TYPES <- c('skene_InN', 'skene_MSN', 'CGE_1', 'CGE_2', 'LGE_1',
                     'LGE_2', 'LGE_4', 'MGE_2', 'MGE_3')

# MAGMA - prepare df
cat('\nPreparing MAGMA data ... \n')
for (LEVEL in c('1', '2')) {

  for (DISORDER in GWAS) {
    
    for (WINDOW in GENE_WINDOW) {
      
      for (DISORDER in GWAS) {
        
        if (LEVEL == '1') {
          
          BF_CORR <- 0.05/6
          
        } else {
          
          BF_CORR <- 0.05/30
          
        }

      MAGMA_DF <- read.table(paste0(MAGMA_DATA_DIR, 'snRNAseq_GE_', DISORDER, '.shi_bc.lvl_', 
                                          LEVEL, '.magma.', WINDOW, '.gsa.out'), header = FALSE) %>%
        janitor::row_to_names(row_number = 1) %>% 
        mutate(VARIABLE = gsub('\\.', '-', VARIABLE)) %>%
        mutate(MAGMA = -log10(as.numeric(P))) %>%
        select(VARIABLE, MAGMA) %>%
        dplyr::rename(Category = VARIABLE) %>%     
        mutate(across('Category', str_replace, 'GE', 'GE-N')) %>%
        mutate(across('Category', str_replace, 'Early_InN', 'IPC')) %>%
        mutate(across('Category', str_replace, '_', '-')) 
      
      MAGMA_PLOT <- ggplot(data = MAGMA_DF, aes(x = MAGMA, y = factor(Category, rev(levels(factor(Category)))), 
                                                fill = '#F8766D')) +
        geom_bar(stat = "identity", color = 'black', position = "dodge") +
        geom_vline(xintercept=-log10(BF_CORR), linetype = "dashed", color = "black") +
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
              axis.text.y  = element_text(colour = "#000000", size = 13)) +
        xlab(expression(-log[10](P))) +
        ylab('Cell type') +
        xlim(0, 8) 
               
      assign(paste0('magma_', DISORDER, '_lvl_', LEVEL, '_', WINDOW, '_df'), MAGMA_DF, envir = .GlobalEnv) 
      assign(paste0('magma_', DISORDER, '_lvl_', LEVEL, '_', WINDOW, '_plot'), MAGMA_PLOT, envir = .GlobalEnv) 
                    
      }
 
    }
  
  }
  
}

# LDSR - prepare df
cat('\nPreparing LDSR data ... \n')
for (DISORDER in GWAS) {
  
  for (LEVEL in c('1', '2')) {
    
    for (WINDOW in GENE_WINDOW) {
      
      if (LEVEL == 1) {
        
        BF_CORR <- 0.05/6
        
        # LDSR                                         
        LDSR_FULL_DF <- read_tsv(paste0(LDSC_DATA_DIR, 'snRNAseq_LDSR_', DISORDER, '_baseline.v1.2_summary.tsv')) %>%
          mutate(LDSR = if_else(`Coefficient_z-score` > 0, -log10(pnorm(`Coefficient_z-score`, lower.tail = FALSE)), 0)) %>%
          separate(Category, into = c('Category', 'Window'), sep = '\\.') %>%
          filter(Category %in% c('LGE', 'MGE', 'CGE', 'Early_InN', 'Progenitor', 'Microglia')) %>%
          mutate(across('Category', str_replace, 'GE', 'GE-N')) %>%
          mutate(across('Category', str_replace, 'Early_InN', 'IPC')) %>%
          mutate(across('Category', str_replace, '_', '-')) %>%
          filter(Window == (!!WINDOW)) 
        
      } else {
        
        BF_CORR <- 0.05/30
        
        LDSR_FULL_DF <- read_tsv(paste0(LDSC_DATA_DIR, 'snRNAseq_LDSR_', DISORDER, '_baseline.v1.2_summary.tsv')) %>%
          mutate(LDSR = if_else(`Coefficient_z-score` > 0, -log10(pnorm(`Coefficient_z-score`, lower.tail = FALSE)), 0)) %>%
          separate(Category, into=c('Category', 'Window'), sep = '\\.') %>%
          filter(str_detect(Category, "_|Other")) %>%
          filter(!str_detect(Category, "Early_InN")) %>%
          filter(Window == (!!WINDOW)) %>%
         mutate(across('Category', str_replace, 'GE', 'GE-N')) %>%
          mutate(across('Category', str_replace, 'Early_InN', 'IPC')) %>%
          mutate(across('Category', str_replace, '_', '-')) 
        
      }
        
      LDSR_DF <- LDSR_FULL_DF %>%
        select(Category, LDSR)
      
      LDSR_PLOT <- ggplot(data = LDSR_DF, aes(x = LDSR, y = factor(Category, rev(levels(factor(Category)))), 
                                                fill = '#F8766D')) +
        geom_bar(stat = "identity", color = 'black', position = "dodge") +
        geom_vline(xintercept=-log10(BF_CORR), linetype = "dashed", color = "black") +
        geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
        theme_bw() +
        ggtitle(paste0(DISORDER, '_', WINDOW)) +
        theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.border = element_rect(colour = "black", size = 1),
              plot.title = element_text(hjust = 0.5, face = 'bold'),
              axis.title.x = element_text(colour = "#000000", size = 14),
              axis.title.y = element_text(colour = "#000000", size = 14),
              axis.text.x  = element_text(colour = "#000000", size = 13, vjust = 0.5),
              axis.text.y  = element_text(colour = "#000000", size = 13)) +
        xlab(expression(-log[10](P))) +
        ylab('Cell type') +
        xlim(0, 8) 
        
      
      assign(paste0('ldsr_', DISORDER, '_lvl_', LEVEL, '_', WINDOW, '_df'), LDSR_DF, envir = .GlobalEnv) 
      assign(paste0('ldsr_', DISORDER, '_lvl_', LEVEL, '_', WINDOW, '_full_df'), LDSR_FULL_DF, envir = .GlobalEnv)
      assign(paste0('ldsr_', DISORDER, '_lvl_', LEVEL, '_', WINDOW, '_plot'), LDSR_PLOT, envir = .GlobalEnv) 
      
    }
    
  }
    
}

# Plot MAGMA and LDSR barplot 
cat('\nCreate plots ... \n')
for (LEVEL in c('1', '2')) {
  
  for (DISORDER in GWAS) {
    
    if (LEVEL == '1') {
      
      BF_CORR <- 0.05/6
      
    } else {
      
      BF_CORR <- 0.05/30
      
    }
    
    PLOT_DF <- left_join(get(paste0('magma_', DISORDER, '_lvl_', LEVEL, '_35UP_10DOWN_df')), 
                         get(paste0('ldsr_', DISORDER, '_lvl_', LEVEL, '_100UP_100DOWN_df')),
                         by = 'Category') %>%
               reshape2::melt()
        
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
            legend.text = element_text(size = 13),
            legend.title = element_blank()) +
      xlab(expression(-log[10](P))) +
      ylab('Cell type') +
      xlim(0, 8) 
      
    assign(paste0(DISORDER, '_magma_ldsr_lvl_', LEVEL, '_plot'), MAGMA_LDSR_PLOT, envir = .GlobalEnv) 
    
  }

}


# MAGMA conditional - prepare df
cat('\nPreparing MAGMA conditional data ... \n')
for (CELL_TYPE in COND_CELL_TYPES) {
  
  for (WINDOW in GENE_WINDOW) { 
  
  MAGMA_DF <- read.table(paste0(MAGMA_COND_DATA_DIR, 'magma_all_sig_and_skene_condition_', 
                                CELL_TYPE, '.', WINDOW, '.gsa.out'), header = TRUE) %>%
    mutate(MAGMA = -log10(as.numeric(P))) %>%
    filter(!grepl(CELL_TYPE, VARIABLE)) %>%
    filter(!grepl('skene', VARIABLE)) %>%
    select(VARIABLE, MAGMA) %>%
    dplyr::rename(Category = VARIABLE) %>%
    dplyr::mutate(across(Category, str_replace, 'GE_', 'GE-N-')) 
  
  DENOM <- if (grepl('skene', CELL_TYPE)) 7 else 6
  
  # if (grepl('skene', CELL_TYPE)) {
  # 
  #   MAGMA_DF <- MAGMA_DF %>%
  #   filter(!grepl('Adult', Category))
  # 
  #   print(MAGMA_DF)
  # 
  #   # } else {
  #   #
  #   #   MAGMA_DF <- MAGMA_DF %>%
  #   #   filter(!grepl(paste0('skene|', CELL_TYPE), Category))
  #   #
  #   }
  
  MAGMA_PLOT <- ggplot(data = MAGMA_DF, aes(x = MAGMA, y = factor(Category, rev(levels(factor(Category)))), 
                                            fill = '#F8766D')) +
    geom_bar(stat = "identity", color = 'black', position = "dodge") +
    geom_vline(xintercept=-log10(0.05/DENOM), linetype = "dashed", color = "black") +
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
    xlim(0, 8) 
  
  
  
  
  assign(paste0('magma_cond_GE_', CELL_TYPE, '_', WINDOW, '_df'), MAGMA_DF, envir = .GlobalEnv) 
  assign(paste0('magma_cond_GE_', CELL_TYPE, '_', WINDOW, '_plot'), MAGMA_PLOT, envir = .GlobalEnv)
  
  }
  
}

# MAGMA downsampled
# MAGMA - prepare df
cat('\nPreparing MAGMA data ... \n')
for (LEVEL in c('1', '2')) {
  
  for (DISORDER in GWAS) {
    
    for (WINDOW in GENE_WINDOW) {
    
      MAGMA_DF <- read.table(paste0(MAGMA_DATA_DIR, 'snRNAseq_GE_', DISORDER, 
                                        '.shi_bc_dwnSmpl.lvl_', LEVEL, '.magma.', 
                                        WINDOW, '.gsa.out'), 
                                 header = FALSE) %>%
      janitor::row_to_names(row_number = 1) %>% 
      mutate(VARIABLE = gsub('\\.', '-', VARIABLE)) %>%
      mutate(MAGMA = -log10(as.numeric(P))) %>%
      select(VARIABLE, MAGMA) %>%
      dplyr::rename(Category = VARIABLE) %>%
      dplyr::mutate(across(Category, str_replace, 'Early_InN', 'IPC')) %>%  
      dplyr::mutate(across(Category, str_replace, '_', '-')) 
    
      MAGMA_PLOT <- ggplot(data = MAGMA_DF, aes(x = MAGMA, y = factor(Category, rev(levels(factor(Category)))), 
                                                   fill = '#F8766D')) +
      geom_bar(stat = "identity", color = 'black', position = "dodge") +
      geom_vline(xintercept=-log10(0.05/6), linetype = "dashed", color = "black") +
      geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
      theme_bw() +
      ggtitle(paste0(DISORDER, "_", WINDOW)) +
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
      xlim(0, 8) 
      
      assign(paste0('magma_dwnSmpl_', DISORDER, '_lvl_', LEVEL, '_', WINDOW, '_df'), MAGMA_DF, envir = .GlobalEnv)   
      assign(paste0('magma_dwnSmpl_', DISORDER, '_lvl_', LEVEL, '_', WINDOW, '_plot'), MAGMA_PLOT, envir = .GlobalEnv) 
    

    }
    
  }
  
}

# LDSR downsampled
# LDSR - prepare df
cat('\nPreparing LDSR data ... \n')
for (DISORDER in GWAS) {
  
  for (LEVEL in c('1', '2')) {
    
    for (WINDOW in GENE_WINDOW) {
      
      if (LEVEL == 1) {
        
        BF_CORR <- 0.05/6
        
        # LDSR                                         
        LDSR_FULL_DF <- read_tsv(paste0(LDSC_DATA_DIR, '/shi_bc_dwnSmpl/snRNAseq_LDSR_', DISORDER, '_baseline.v1.2_summary.tsv')) %>%
          mutate(LDSR = if_else(`Coefficient_z-score` > 0, -log10(pnorm(`Coefficient_z-score`, lower.tail = FALSE)), 0)) %>%
          separate(Category, into = c('Category', 'Window'), sep = '\\.') %>%
          filter(Category %in% c('LGE', 'MGE', 'CGE', 'Early_InN', 'Progenitor', 'Microglia')) %>%
          dplyr::mutate(across(Category, str_replace, 'Early_InN', 'IPC')) %>%
          dplyr::mutate(across(Category, str_replace, 'GE', 'GE-N')) %>%
          filter(Window == (!!WINDOW)) 
        
      } else {
        
        BF_CORR <- 0.05/30
        
        LDSR_FULL_DF <- read_tsv(paste0(LDSC_DATA_DIR, '/shi_bc_dwnSmpl/snRNAseq_LDSR_', DISORDER, '_baseline.v1.2_summary.tsv')) %>%
          mutate(LDSR = if_else(`Coefficient_z-score` > 0, -log10(pnorm(`Coefficient_z-score`, lower.tail = FALSE)), 0)) %>%
          separate(Category, into=c('Category', 'Window'), sep = '\\.') %>%
          filter(str_detect(Category, "_|Other")) %>%
          filter(!str_detect(Category, "Early_InN")) %>%
          dplyr::mutate(across(Category, str_replace, 'GE', 'GE-N')) %>%
          dplyr::mutate(across(Category, str_replace, '_', '-')) %>%
          filter(Window == (!!WINDOW)) 
        
      }
      
      LDSR_DF <- LDSR_FULL_DF %>%
        select(Category, LDSR)
      
      LDSR_PLOT <- ggplot(data = LDSR_DF, aes(x = LDSR, y = factor(Category, rev(levels(factor(Category)))), 
                                              fill = '#F8766D')) +
        geom_bar(stat = "identity", color = 'black', position = "dodge") +
        geom_vline(xintercept=-log10(BF_CORR), linetype = "dashed", color = "black") +
        geom_vline(xintercept=-log10(0.05), linetype = "dotted", color = "black") +
        theme_bw() +
        ggtitle(paste0(DISORDER, '_', WINDOW)) +
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
        xlim(0, 8) 
      
      
      assign(paste0('ldsr_dwnSmpl_', DISORDER, '_lvl_', LEVEL, '_', WINDOW, '_df'), LDSR_DF, envir = .GlobalEnv) 
      assign(paste0('ldsr_dwnSmpl_', DISORDER, '_lvl_', LEVEL, '_', WINDOW, '_full_df'), LDSR_FULL_DF, envir = .GlobalEnv)
      assign(paste0('ldsr_dwnSmpl_', DISORDER, '_lvl_', LEVEL, '_', WINDOW, '_plot'), LDSR_PLOT, envir = .GlobalEnv) 
      
    }
    
  }
  
}

# Produce final plots for paper. ------------------------------------------------------
legend <- get_legend(SCZ_magma_ldsr_lvl_1_plot)

# Level 1 - Gene set enrichment barplots
figure_2 <- plot_grid(SCZ_magma_ldsr_lvl_1_plot + NoLegend(), 
                      HEIGHT_magma_ldsr_lvl_1_plot + NoLegend(), 
                      legend, ncol = 3, rel_widths = c(1, 1, 0.5),  
                      labels = c('A', 'B', ''), label_size = 20)

# Level 2 - Gene set enrichment barplots
figure_4 <- plot_grid(SCZ_magma_ldsr_lvl_2_plot + NoLegend(), 
                      HEIGHT_magma_ldsr_lvl_2_plot+ NoLegend(), 
                      legend, ncol = 3, rel_widths = c(1, 1, 0.5), 
                      labels = c('A', 'B', ''), label_size = 20)

# Figure S1 - Gene window plots
gene_windows_lvl1_plot <- plot_grid(magma_SCZ_lvl_1_10UP_10DOWN_plot + NoLegend() + ggtitle(paste0('MAGMA 10UP-10DOWN')), 
                                    magma_SCZ_lvl_1_35UP_10DOWN_plot + NoLegend() + ggtitle(paste0('MAGMA 35UP-10DOWN')), 
                                    magma_SCZ_lvl_1_100UP_100DOWN_plot + NoLegend() + ggtitle(paste0('MAGMA 100UP-100DOWN')),
                                    ldsr_SCZ_lvl_1_10UP_10DOWN_plot + NoLegend() + ggtitle(paste0('SLDSR 10UP-10DOWN')), 
                                    ldsr_SCZ_lvl_1_35UP_10DOWN_plot + NoLegend() + ggtitle(paste0('SLDSR 35UP-10DOWN')), 
                                    ldsr_SCZ_lvl_1_100UP_100DOWN_plot + NoLegend() + ggtitle(paste0('SLDSR 100UP-100DOWN')), 
                                    ncol = 3, labels = c('AUTO'), label_size = 20)

# Figure S2 - downsampled plots
dwnSmple_plot <- plot_grid(magma_SCZ_lvl_1_35UP_10DOWN_plot + NoLegend() + ggtitle(paste0('MAGMA Level 1')), 
                           ldsr_dwnSmpl_SCZ_lvl_1_100UP_100DOWN_plot + NoLegend() + ggtitle(paste0('SLDSR Level 1')),
                           magma_dwnSmpl_SCZ_lvl_2_35UP_10DOWN_plot + NoLegend() + ggtitle(paste0('MAGMA Level 2')), 
                           ldsr_dwnSmpl_SCZ_lvl_2_100UP_100DOWN_plot + NoLegend() + ggtitle(paste0('SLDSR Level 2')),
                           ncol = 2, labels = c('AUTO'), label_size = 20, rel_heights = c(1, 3))

# Figure S7 - Within Level 2 cell conditional
plot_grid(magma_cond_GE_CGE_1_35UP_10DOWN_plot + NoLegend() + ggtitle(paste0('CGE-N-1')), 
          magma_cond_GE_CGE_2_35UP_10DOWN_plot + NoLegend() + ggtitle(paste0('CGE-N-2')),
          magma_cond_GE_LGE_1_35UP_10DOWN_plot + NoLegend() + ggtitle(paste0('LGE-N-1')), 
          magma_cond_GE_LGE_2_35UP_10DOWN_plot + NoLegend() + ggtitle(paste0('LGE-N-2')),
          magma_cond_GE_LGE_4_35UP_10DOWN_plot + NoLegend() + ggtitle(paste0('LGE-N-4')),
          magma_cond_GE_MGE_2_35UP_10DOWN_plot + NoLegend() + ggtitle(paste0('MGE-N-2')),
          magma_cond_GE_MGE_3_35UP_10DOWN_plot + NoLegend() + ggtitle(paste0('MGE-N-3')),
          ncol = 3, labels = c('AUTO'), label_size = 20)

# Figure S8 - Adult conditional
plot_grid(magma_cond_GE_skene_InN_35UP_10DOWN_plot + NoLegend() + ggtitle(paste0('Adult InN')), 
          magma_cond_GE_skene_MSN_35UP_10DOWN_plot + NoLegend() + ggtitle(paste0('Adult MSN')),
          ncol = 2, labels = c('AUTO'), label_size = 20)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------


