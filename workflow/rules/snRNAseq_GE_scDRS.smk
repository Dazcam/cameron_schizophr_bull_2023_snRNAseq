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
rule scDRS_create_zscore_files:
    input:   "../results/magma/snRNAseq_GE_SCZ.magma.genes.out"
    output:  "../results/scDRS/scDRS_genewise_Z.tsv"
    log:     "../results/logs/scDRS/snRNAseq.GE.create_zscore_files.log"
    shell:
             """
             
             awk -v OFS="\t" '$1=$1' {input} | cut -f 1,8 | sed 's/ZSTAT/SCZ/g' > {output} &> {log}

             """

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

rule create_covariate_file:
    # Works manually but not in snakerule yet???
    input:   "../results/h5ad_objects/shi2021_filt.h5ad"
    output:  "../results/scDRS/scDRS_covariates.tsv"
    log:     "../results/logs/magma/snRNAseq.GE.create_covariates.log"
    script:   "../scripts/snRNAseq_GE_create_covariates.py"         


rule scDRS_compute_score:
    input:   rna_data = "../results/h5ad_objects/shi2021_filt.h5ad", 
             gs = "../results/scDRS/scDRS_genewise_Z_top1K.gs",
             cov = "../results/scDRS/scDRS_covariates.tsv"
    output:  "../results/scDRS/scDRS_SCZ.score.gz"
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
    input:   rna_data =	"../results/h5ad_objects",
             gwas_scores = "../results/scDRS/SCZ.full_score.gz"
    output:  "../results/scDRS/SCZ.scdrs_group.level1class"
    params:  "../results/scDRS/"
    log:     "../results/logs/magma/snRNAseq.GE.group_level_stats.log"
    shell:
             """

             scdrs perform-downstream \
               --h5ad-file {input.rna_data} \
               --score-file {input.gwas_scores} \
               --out-folder {params} \
               --group-analysis level1class \
               --flag-filter-data True \
               --flag-raw-count True &> {log}
      
       	      """

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
