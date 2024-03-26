# Genetic implication of prenatal GABAergic and cholinergic neuron development in susceptibility to schizophrenia (2024)

This project was carried out in the Division of Psychological Medicine and Clinical Neurosciences (DPMCN). The paper is [here]() ADD WHEN PUBLISHED. 

***

## **snRNAseq Data**

snRNAseq data for this study this study is taken from [Shi et al. (2021)](https://www.science.org/doi/10.1126/science.abj6641):

+ [GEO - GSE135827](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135827)
+ [Metadata](https://www.science.org/doi/suppl/10.1126/science.abj6641/suppl_file/science.abj6641_tables_s2_to_s9.zip)
+ [GeX Matrix](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135827/suppl/GSE135827%5FGE%5Fmat%5Fraw%5Fcount%5Fwith%5Fweek%5Finfo%2Etxt%2Egz)

Details on the snATACseq data  / analyses for this study can be found [here](https://github.com/Dazcam/cameron_schizophr_bull_2023_snATACseq).

***

## **GWAS Data**

See the following papers for GWAS data access:

+ [Schizophrenia](https://figshare.com/ndownloader/files/28169757)
+ [Autism](https://figshare.com/ndownloader/files/28169292)
+ [Major Depressive Disorder]() - Permission required at time of access
+ [ADHD](https://figshare.com/ndownloader/files/40036684)
+ [Bipolar Disorder]() - Permission required at time of access
+ [Height](https://portals.broadinstitute.org/collaboration/giant/images/6/63/Meta-analysis_Wood_et_al%2BUKBiobank_2018.txt.gz)

***

**Scripts**

1. [snRNAseq_GE_get_shi_data.sh](scripts/snRNAseq_GE_get_shi_data.sh) - Download Shi data
2. [snRNAseq_GE_prep_shi_data_for_Seurat.R](scripts/snRNAseq_GE_prep_shi_data_for_Seurat.R) - Clean and prepare Shi data
3. [snRNAseq_GE_shi_seurat.R](scripts/snRNAseq_GE_shi_seurat.R) - Generate Seurat objects (with batch correction) for cluster levels 1 and 2
4. [snRNAseq_GE_prep_enrich_test_files.R](scripts/snRNAseq_GE_prep_enrich_test_files.R) - Generate specificity scores and top 10% genesets for enrichment testing



