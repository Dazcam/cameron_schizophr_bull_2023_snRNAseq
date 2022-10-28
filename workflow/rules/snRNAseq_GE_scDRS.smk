# -------------------------------------------------------------------------------------
#
#
#    scDRS analyses
#
#
# -------------------------------------------------------------------------------------


# ---------  SET SMK PARAMS  ----------
configfile: "../config/config.yaml"

# -------------  RULES  ---------------
rule scDRS_create_h5ad_and_cov_files:
    input:   "../results/h5ad_objects/{H5AD_OBJ}.h5ad"
    output:  "../results/h5ad_objects/{H5AD_OBJ}.qc.h5ad",
             "../results/scDRS/scDRS_{H5AD_OBJ}/scDRS_{H5AD_OBJ}_covariates.tsv"
    log:     "../results/logs/scDRS/snRNAseq.GE.create_h5ad_and_cov_{H5AD_OBJ}.log"
    script:   "../scripts/snRNAseq_GE_prepare_data_for_scDRS.py"


rule scDRS_create_zscore_files:
    input:   "../results/magma/snRNAseq_GE_SCZ.magma.genes.out", 
             "../resources/refs/NCBI37.3.gene.loc.extendedMHCexcluded.txt"
    output:  "../results/scDRS/scDRS_genewise_Z.tsv"
    log:     "../results/logs/scDRS/snRNAseq.GE.create_zscore_files.log"
    script:  "../scripts/snRNAseq_GE_create_z_scores.py"

rule scDRS_munge_gs:
    input:   "../results/scDRS/scDRS_genewise_Z.tsv"
    output:  "../results/scDRS/scDRS_genewise_Z_top1K.gs"
    log:     "../results/logs/scDRS/snRNAseq.GE.munge_gs.log"
    shell:
             """
             
             scdrs munge-gs \
                   --out-file {output} \
                   --zscore-file {input} \
                   --weight zscore \
                   --n-max 1000 &> {log}

             """        

rule scDRS_compute_score:
    input:   rna_data = "../results/h5ad_objects/{H5AD_OBJ}.qc.h5ad", 
             gs = "../results/scDRS/scDRS_genewise_Z_top1K.gs",
             cov = "../results/scDRS_{H5AD_OBJ}/scDRS_{H5AD_OBJ}_covariates.tsv"
    output:  "../results/scDRS/scDRS_{H5AD_OBJ}/SCZ.full_score.gz"
    params:  "../results/scDRS/scDRS_{H5AD_OBJ}/"
    log:     "../results/logs/scDRS/snRNAseq.GE.compute.scores_{H5AD_OBJ}.log"
    shell:
             """
          
             scdrs compute-score \
               --h5ad-file {input.rna_data} \
               --h5ad-species human \
               --gs-file {input.gs} \
               --gs-species human \
               --cov-file {input.cov} \
               --flag-filter-data True \
               --flag-raw-count True \
               --flag-return-ctrl-raw-score False \
               --flag-return-ctrl-norm-score True \
               --out-folder {params} &> {log}
      
              """

rule scDRS_group_level_stats_all:
    input:   rna_data =	"../results/h5ad_objects/shi_bc.qc.h5ad",
             gwas_scores = "../results/scDRS/scDRS_shi_bc/SCZ.full_score.gz"
    output:  "../results/scDRS/scDRS_shi_bc/SCZ.scdrs_group.cluster_level_1"
    params:  "../results/scDRS/scDRS_shi_bc/"
    log:     "../results/logs/scDRS/snRNAseq.GE.group_level_stats_shi_bc.log"
    shell:
             """

             scdrs perform-downstream \
               --h5ad-file {input.rna_data} \
               --score-file {input.gwas_scores} \
               --out-folder {params} \
               --group-analysis cluster_level_1 \
               --flag-filter-data True \
               --flag-raw-count True &> {log}
      
       	      """

rule scDRS_group_level_stats_subclust:
    input:   rna_data = "../results/h5ad_objects/{H5AD_OBJ_SUBCLUST}.qc.h5ad",
             gwas_scores = "../results/scDRS/scDRS_{H5AD_OBJ_SUBCLUST}/SCZ.full_score.gz"
    output:  "../results/scDRS/scDRS_{H5AD_OBJ_SUBCLUST}/SCZ.scdrs_group.cluster_level_2"
    params:  "../results/scDRS/scDRS_{H5AD_OBJ_SUBCLUST}/"
    log:     "../results/logs/scDRS/snRNAseq.GE.group_level_stats_{H5AD_OBJ_SUBCLUST}.log"
    shell:
             """

             scdrs perform-downstream \
               --h5ad-file {input.rna_data} \
               --score-file {input.gwas_scores} \
               --out-folder {params} \
               --group-analysis cluster_level_2 \
               --flag-filter-data True \
               --flag-raw-count True &> {log}

              """

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
