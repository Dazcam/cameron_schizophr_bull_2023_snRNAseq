#--------------------------------------------------------------------------------------
#
#    snRNAseq GE - plots for paper
#
#--------------------------------------------------------------------------------------

## Load packages  ---------------------------------------------------------------------
if (!require("Require")) install.packages("Require")
Require::Require(c("tidyverse", "readxl", "ComplexHeatmap", "pheatmap", "ArchR",
                   "Seurat", "gtools", "plyr", "cowplot", "scCustomize"))

## Load env variables  ----------------------------------------------------------------
SCRIPT_RNA_DIR <- '~/Desktop/fetal_brain_snRNAseq_GE_270922/workflow/scripts/'
SCRIPT_ATAC_DIR <- '~/Desktop/fetal_brain_snATACseq_V3_010323/workflow/scripts/'
RESULTS_RNA_DIR <- '~/Desktop/fetal_brain_snRNAseq_GE_270922/results/'
RESULTS_ATAC_DIR <- '~/Desktop/fetal_brain_snATACseq_V3_010323/results/'
ARCHR_DIR <- paste0(RESULTS_ATAC_DIR, '02ARCHR/')
FRAGS_DIR <- paste0(RESULTS_ATAC_DIR, '04FRAGMENT_FILES/')
PEAKS_DIR <- paste0(RESULTS_ATAC_DIR, '05PEAKS/')
R_DIR <- paste0(RESULTS_RNA_DIR, '01R_objects/')

# ATAC gene lists
MGE_GENES <- c('CRABP1', 'NKX2-1', 'ETV1', 'NFIA', # Striatal and CRT Ns
               'MEF2C', 'MAF', 'ARX', 'ZEB2', # Cortical interneurons
               'ZFHX3', # Subpallial and GABA Ns also tested 'NR2F1', 'NR2F2', 
               'LHX8', 'CNTNAP2', 'GRIA2', 'ZIC1', # Subpallial Cholin Ns also tested 'ISL1', 'GBX2',
               'MDK', 'LHX6', 'SOX6', 'CXCR4', # Rest from Shi fig 5B
               'ERBB4', 'ANGPT2', 
               'SST'

)

LGE_GENES <- c(#'MEIS2', 'TLE4', 'FOXP1', # 'ZNF503' (duplicate), # Striatal potential
  'CHD7', 'DLX5', 'ID2', 'PAX6', # OB precursor
  'PENK', 'CXCL12', 'SP9', # D2 MSN precursor
  'TAC1', 'EBF1', 'ISL1', 'ZNF503', 'PDYN', # D1 MSN precursor (PDYN)
  'TSHZ1', 'FOXP2', 'PBX3', 'ERBB4' #, # D1 MSN precursor (TSHZ1)
  #'SIX3', 'ID4', 'ZFHX4', 'PHLDA1', # Rest from Shi fig 4H
  #'GAP43', 'ZNF521', 'NTRK3',
  #'LHFP', 'ZADH2', 'PBX1',
  #'ADRA2A', 'DRD1', 'DRD2'
)

# LGE_GENES <- c('FOXP1', 'FOXP2', 'ISL1', 'SERTAD4', 'ZNF503',
#                'SIX3', 'ZFHX3', 'SCGN', 'PBX1', 'ZNF521',
#                'PDYN', 'EBF1', 'GRIA2', 'CNTNAP2', 'ZFHX3')
CGE_GENES <- c(
  'NRIP3', 'IER2', 'PRKCA', 'ERBB4',  #CGE-InN-1
  'VSTM2A', 'BEX2', 'PRKCA', 'ERBB4', 'ARL6IP5', # CGE-InN-2
  'ARL4D', 'HIGD2A', 'GNG5', 'SP9', 'STMN4','PBX1', 'MDK', 'SOX6', 'HES6',
  'GAD2', 'CALB2' #CGE-InN-3
)

PROGENITOR_GENES <- c('SOX6', 'DLX6', 'SFTA3', 'NKX2-1', # 'LHX6', 
                      'RIC3', 'MBIP', 'OLIG2', 'PCDH17', # 'GSX1',
                      'NR2F1', 'NR2F2', 'NTRK2', 'DLX5', 'RBP1', 
                      'ZFHX3', 'MEIS2', # 'VIP', 'CALB2', 'SCGN',
                      'SIX3', 'PHLDA1', 'GLI3', 'LGALS1', # 'ZNF503',
                      'PAX6', 'GSX2', 'HEY1'
                      
)

IPC_GENES <- c('DLL1', 'DLL3', 'CCND2', 'GAD1', 'GAD2') # 'KI67', 'EOMES'

LEVEL_1_MARKERS <- c("CALB2", "SCGN", "PCDH9", "ANKS1B",   # CGE
                     "FOXP1", "ZNF503", "SERTAD4", "ISL1", # LGE
                     "LHX6", "NXPH1",                      # MGE
                     "GAD2",                               # InN
                     "HES1", 'OLIG2', "PAX6", "GSX2", 'PTPRZ1', # Progenitor
                     "DLL1", "DLL3", "CCND2",              # IPC
                     "SPI1", "CD68")                       # MG



# RNA gene lists
greens <- c('#3CBB75FF', '#00FF00A5','#006400', '#B7FFB7', '#10A53DFF',
            '#95D840FF', '#9DC183',  '#708238', '#55C667FF', '#73D055FF',
            '#567D46')

ExN_blues <- c('#76B5C5', '#00BDD2', '#CEE5FD', '#00B6EB', '#ABDBE3',
              '#1E81B0', '#B8D2EB', '#779CBA')

reds <- c('#FAA0A0', '#FF5959', '#F75151', '#EF0029', '#D2042D')

reds <- c('#B200ED',  '#DCBEFF', '#6F2DA8')





## Load functions  --------------------------------------------------------------------
source(paste0(SCRIPT_ATAC_DIR, 'snATACseq_functions.R'))
source(paste0(SCRIPT_RNA_DIR, 'snRNAseq_GE_plots_magma_and_ldsr_plots.R'))
source(paste0(SCRIPT_ATAC_DIR, 'snATACseq_plots_LDSR_barplots.R'))

## RNA seq ----------------------------------------------------------------------------
## Load Seurat Objects
seurat.shi.bc <- readRDS(paste0(R_DIR, 'seurat_shi_bc.rds'))
seurat.shi.bc_LGE <- readRDS(paste0(R_DIR, 'seurat_shi_bc_LGE.rds'))
seurat.shi.bc_MGE <- readRDS(paste0(R_DIR, 'seurat_shi_bc_MGE.rds'))
seurat.shi.bc_CGE <- readRDS(paste0(R_DIR, 'seurat_shi_bc_CGE.rds'))
seurat.shi.bc_Progenitor <- readRDS(paste0(R_DIR, 'seurat_shi_bc_Progenitor.rds'))
seurat.shi.bc_Early_InN <- readRDS(paste0(R_DIR, 'seurat_shi_bc_Early_InN.rds'))
seurat.shi.bc_dwnSmpl_lvl1 <- readRDS(paste0(R_DIR, 'seurat_shi_bc_dwnSmpl_lvl1.rds'))
seurat.shi.bc_dwnSmpl_lvl2 <- readRDS(paste0(R_DIR, 'seurat_shi_bc_dwnSmpl_lvl2.rds'))

## Rename clusters
for (SUFFIX in c('', '_LGE', '_MGE', '_CGE', '_Progenitor',
                 '_Early_InN', '_dwnSmpl_lvl1', '_dwnSmpl_lvl2' )) {

  cat('\nRunning seurat.shi.bc', SUFFIX, sep = '')

  SEURAT_OBJ <- get(paste0('seurat.shi.bc', SUFFIX))

  if (SUFFIX == '' |  SUFFIX == '_dwnSmpl_lvl1' | SUFFIX == '_dwnSmpl_lvl2') {

    CLUST = 'cluster_level_1'
    NEW_NAMES <- SEURAT_OBJ[[]] %>%
      as_tibble() %>%
      dplyr::select(cluster_level_1) %>%
      dplyr::mutate(across(CLUST, str_replace, 'Early_InN', 'IPC')) %>%
      dplyr::mutate(across(CLUST, str_replace, 'GE', 'GE-N')) %>%
      pull(cluster_level_1)

  } else {

    CLUST = 'cluster_level_2'
    NEW_NAMES <-SEURAT_OBJ[[]] %>%
      as_tibble() %>%
      dplyr::select(cluster_level_2) %>%
      dplyr::mutate(across(CLUST, str_replace, 'Early_InN', 'IPC')) %>%
      dplyr::mutate(across(CLUST, str_replace, 'GE', 'GE-N')) %>%
      dplyr::mutate(across(CLUST, str_replace, '_', '-')) %>%
      pull(cluster_level_2)


    }


  #

  SEURAT_OBJ[[paste0(CLUST, '_new')]] <- NEW_NAMES
  print(table(SEURAT_OBJ[[paste0(CLUST, '_new')]]))
  assign(paste0('seurat.shi.bc', SUFFIX), SEURAT_OBJ, .GlobalEnv)

}

## Cell count table
cnts <- data.frame(Region_cell_type  = c("Ganglionic eminence", "CGE-N", "LGE-N", "MGE-N", "Progenitor", "IPC"),
                   Level = c(1, 2, 2, 2, 2, 2),
                   Cell_cnt = c(ncol(seurat.shi.bc), ncol(seurat.shi.bc_CGE), ncol(seurat.shi.bc_LGE), 
                                ncol(seurat.shi.bc_MGE), ncol(seurat.shi.bc_Progenitor), ncol(seurat.shi.bc_Early_InN)),
                   Gene_cnt = c(nrow(seurat.shi.bc), nrow(seurat.shi.bc_CGE), nrow(seurat.shi.bc_LGE), 
                                nrow(seurat.shi.bc_MGE), nrow(seurat.shi.bc_Progenitor), nrow(seurat.shi.bc_Early_InN))
)


## FINAL PLOTS FOR RNA SEQ PAPER (EXCLUDING BAR PLOTS)  ---------------------
# Fig 1 - L1 cluster UMAP and Vln  -----
seurat.shi.bc$cluster_level_1_new <- factor(seurat.shi.bc$cluster_level_1_new,
                                            level=c('LGE-N', 'MGE-N', 'CGE-N',
                                                    'IPC', 'Microglia', 'Progenitor'))
figure_1A <- DimPlot_scCustom(seurat.shi.bc, group.by = 'cluster_level_1_new',
                              DiscretePalette_scCustomize(num_colors = 26, palette = "ditto_seq"),
                              label = TRUE, repel = TRUE) +
  theme(legend.position = "none") +
  ggtitle("") 

figure_1B <- VlnPlot(seurat.shi.bc, LEVEL_1_MARKERS, stack = TRUE, flip = TRUE,
                    cols = DiscretePalette_scCustomize(num_colors = 26, palette = "ditto_seq"),
                    same.y.lims = TRUE, fill.by = 'ident', group.by = 'cluster_level_1_new') +
  theme(axis.title.x=element_blank(), legend.position = "none")

plot_grid(figure_1A, figure_1B, labels = 'AUTO', label_size = 20)

# Fig S5-S9 - L2 cluster UMAPs and Vlns -----
for (SUFFIX in c('_CGE', '_LGE', '_MGE', '_Progenitor', '_Early_InN')) {

  cat('\nRunning seurat.shi.bc', SUFFIX, sep = '')
  SEURAT_OBJ <- get(paste0('seurat.shi.bc', SUFFIX))

  GET_MARKERS <- function(x) { # Instead of nested ifelse
    switch(x,
           "_CGE" = CGE_GENES,
           "_LGE" = LGE_GENES,
           "_MGE" = MGE_GENES,
           "_Progenitor" = PROGENITOR_GENES,
           "_Early_InN" = IPC_GENES,
           stop("Invalid SUFFIX value")
    )

  }

  MARKERS <- GET_MARKERS(SUFFIX)

  UMAP_PLOT <- DimPlot_scCustom(SEURAT_OBJ, group.by = 'cluster_level_2_new',
                   DiscretePalette_scCustomize(num_colors = 26, palette = "ditto_seq"),
                   label = TRUE) + # Need repel for Prog only otherwise Prog-10 label is cut in image
    theme(legend.position = "none") +
    ggtitle("")

  VLN_PLOT <- VlnPlot(SEURAT_OBJ, MARKERS, stack = TRUE, flip = TRUE,
                       cols = DiscretePalette_scCustomize(num_colors = 26, palette = "ditto_seq"),
                       same.y.lims = TRUE, fill.by = 'ident', group.by = 'cluster_level_2_new') +
    theme(axis.title.x=element_blank(), legend.position = "none")
  
  assign(paste0('cluster', SUFFIX, '_umap'), UMAP_PLOT, .GlobalEnv)
  assign(paste0('cluster', SUFFIX, '_vln'), VLN_PLOT, .GlobalEnv)
  

}

# Fig S5 to S9
fig_S5 <- plot_grid(cluster_LGE_umap, cluster_LGE_vln, labels = 'AUTO', label_size = 20)
fig_S6 <-plot_grid(cluster_MGE_umap, cluster_MGE_vln, labels = 'AUTO', label_size = 20)
fig_S7 <-plot_grid(cluster_CGE_umap, cluster_CGE_vln, labels = 'AUTO', label_size = 20)
fig_S8 <-plot_grid(cluster_Progenitor_umap, cluster_Progenitor_vln, labels = 'AUTO', label_size = 20)
fig_S9 <-plot_grid(cluster_Early_InN_umap, cluster_Early_InN_vln, labels = 'AUTO', label_size = 20)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
