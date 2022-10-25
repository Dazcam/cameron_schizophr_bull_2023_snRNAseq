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
    input:   "../results/h5ad_objects/shi.bc.h5ad"
    output:  "../results/h5ad_objects/shi.bc.qc.h5ad",
             "../results/scDRS/scDRS_covariates.tsv"
    log:     "../results/logs/magma/snRNAseq.GE.create_h5ad_and_cov.log"
    script:   "../scripts/snRNAseq_GE_prepare_data_for_scDRS.py"


rule scDRS_create_zscore_files:
    input:   "../results/magma/snRNAseq_GE_SCZ.magma.genes.out", 
             "../resources/refs/NCBI37.3.gene.loc.extendedMHCexcluded.txt"
    output:  "../results/scDRS/scDRS_genewise_Z.tsv"
    log:     "../results/logs/scDRS/snRNAseq.GE.create_zscore_files.log"
    script:  "../scripts/snRNAseq_GE_create_z_scores.py"

##  Obsolete if rule above works: issue with ENTREZ IDS
#rule scDRS_create_zscore_files:
#    input:   "../results/magma/snRNAseq_GE_SCZ.magma.genes.out"
#    output:  "../results/scDRS/scDRS_genewise_Z.tsv"
#    log:     "../results/logs/scDRS/snRNAseq.GE.create_zscore_files.log"
#    shell:
#             """
             
#             awk -v OFS="\t" '$1=$1' {input} | cut -f 1,8 | sed 's/ZSTAT/SCZ/g' > {output} &> {log}

#             """

rule scDRS_munge_gs:
    input:   "../results/scDRS/scDRS_genewise_Z.tsv"
    output:  "../results/scDRS/scDRS_genewise_Z_top1K.gs"
    log:     "../results/logs/magma/snRNAseq.GE.munge_gs.log"
    shell:
             """
             
             scdrs munge-gs \
                   --out-file {output} \
                   --zscore-file {input} \
                   --weight zscore \
                   --n-max 1000 &> {log}

             """        

##  Obsolete if h5ad rule works: issue with compute_score due to
#rule create_covariate_file:
#    # Works manually but not in snakerule yet???
#    input:   "../results/h5ad_objects/shi2021_filt.h5ad"
#    output:  "../results/scDRS/scDRS_covariates.tsv"
#    log:     "../results/logs/magma/snRNAseq.GE.create_covariates.log"
#    script:   "../scripts/snRNAseq_GE_create_covariates.py"         


rule scDRS_compute_score:
    input:   rna_data = "../results/h5ad_objects/shi.bc.qc.h5ad", 
             gs = "../results/scDRS/scDRS_genewise_Z_top1K.gs",
             cov = "../results/scDRS/scDRS_covariates.tsv"
    output:  "../results/scDRS/SCZ.full_score.gz"
    params:  "../results/scDRS/"
    log:     "../results/logs/magma/snRNAseq.GE.compute.scores.log"
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

rule scDRS_group_level_stats:
    input:   rna_data =	"../results/h5ad_objects/shi.bc.qc.h5ad",
             gwas_scores = "../results/scDRS/SCZ.full_score.gz"
    output:  "../results/scDRS/SCZ.scdrs_group.cluster_level_1"
    params:  "../results/scDRS/"
    log:     "../results/logs/magma/snRNAseq.GE.group_level_stats.log"
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

#rule perform_downstream:
#    input:   rna_data = "../results/h5ad_objects",
#             scores = "../results/scDRS/SCZ.full_score.gz"
#    output:  "../results/scDRS/<trait>.scdrs_group.{CELL_TYPE}"
#    params:  "../results/scDRS/"
#    log:     "../results/logs/magma/snRNAseq.GE.compute.scores.log"
#    shell:
#             """

#             scdrs perform-downstream \
#               --h5ad-file {input.rna_data} \
#               --score-file {input.scores} \
#               --out-folder {params} \
#               --group-analysis {wildcards.CELL_TYPE} \
#               --corr-analysis causal_variable,non_causal_variable,covariate\
#               --gene-analysis \
#               --flag-filter-data True \
#               --flag-raw-count True

#              """



# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
