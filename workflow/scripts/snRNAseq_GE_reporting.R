# -------------------------------------------------------------------------------------
#
#    snRNAseq - Reporting
#
# -------------------------------------------------------------------------------------

##  Set Variables  --------------------------------------------------------------------
SCRIPT_DIR <- '~/Desktop/fetal_brain_snRNAseq_GE_270922/workflow/scripts/'

##  Load objects  --------------------------------------------------------------------
source(paste0(SCRIPT_DIR, 'snRNAseq_GE_scDRS_visualisation.R'))

# Create markdown file
OUT_DIR <- '~/Desktop/fetal_brain_snRNAseq_GE_270922/results/'
MARKDOWN_FILE <- paste0(SCRIPT_DIR, 'snRNAseq_GE_reporting.Rmd')
REPORT_DIR <- paste0(OUT_DIR, 'rmarkdown_reports/')
REPORT_FILE <- 'snRNAseq_GE_reporting.html'


## scDRS
scDRS_all_umap <- plot_grid(scDRS_all_SCZ_group_umap, scDRS_all_HEIGHT_group_umap, ncol = 1)
scDRS_MGE_umap <- plot_grid(scDRS_MGE_SCZ_group_umap, scDRS_MGE_HEIGHT_group_umap, ncol = 1)
scDRS_LGE_umap <- plot_grid(scDRS_LGE_SCZ_group_umap, scDRS_LGE_HEIGHT_group_umap, ncol = 1)
scDRS_CGE_umap <- plot_grid(scDRS_CGE_SCZ_group_umap, scDRS_CGE_HEIGHT_group_umap, ncol = 1)
scDRS_Progenitor_umap <- plot_grid(scDRS_Progenitor_SCZ_group_umap, scDRS_Progenitor_HEIGHT_group_umap, ncol = 1)

render(MARKDOWN_FILE, output_file = REPORT_FILE, output_dir = REPORT_DIR)
