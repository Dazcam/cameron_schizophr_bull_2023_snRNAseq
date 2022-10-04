# -------------------------------------------------------------------------------------
#
#
#    MAGMA analyses
#
#
# -------------------------------------------------------------------------------------


# ---------  SET SMK PARAMS  ----------
configfile: "../config/config.yaml"

# -------------  RULES  ---------------
rule create_zscore_files:
    input:   "../results/magma/snRNAseq_GE_SCZ.magma.genes.out"
    output:  "../results/scDRS/scDRS_genewise_Z.tsv"
    log:     "../results/logs/scDRS/snRNAseq.GE.create_zscore_files.log"
    shell:
             """
             
             awk -v OFS="\t" '$1=$1' test | cut -f 1,8 | sed 's/ZSTAT/SCZ/g' > {output}

             """

rule munge_gs:
    input:   "../results/scDRS/scDRS_genewise_Z.tsv"
    output:  "../results/scDRS/scDRS_genewise_Z_top1K.gs"
    log:     "../results/logs/magma/snRNAseq.GE.munge_gs.log"
    shell:
             """
             
             scdrs munge-gs \
                   --out-file {output} \
                   --zscore-file {input} \
                   --weight zscore \
                   --n-max 1000

             """        

rule compute_score:
    input:   rna_data = "../results/h5ad_objects", 
             gs = "../results/scDRS/scDRS_genewise_Z_top1K.gs"
    output:  "../results/scDRS/scDRS_genewise_Z_top1K.gs"
    params:  "../results/scDRS/"
    log:     "../results/logs/magma/snRNAseq.GE.compute.scores.log"
    shell:
             """
          
             scdrs compute-score \
               --h5ad-file {input.rna_data} \
               --h5ad-species human \
               --gs-file {input.gs} \
               --gs-species human \
               --cov-file data/cov.tsv \
               --flag-filter-data True \
               --flag-raw-count True \
               --flag-return-ctrl-raw-score False \
               --flag-return-ctrl-norm-score True \
               --out-folder {params}
      
              """

rule group_level_stats:
    input:   rna_data =	"../results/h5ad_objects",
             gwas.scores = "../results/scDRS/SCZ.full_score.gz"
    output:  "../results/scDRS/SCZ.scdrs_group.level1class"
    params:  "../results/scDRS/"
    log:     "../results/logs/magma/snRNAseq.GE.group_level_stats.log"
    shell:
             """

             scdrs perform-downstream \
               --h5ad-file {input.rna_data} \
               --score-file {input.gwas.scores} \
               --out-folder {params} \
               --group-analysis level1class \
               --flag-filter-data True \
               --flag-raw-count True	  
      
       	      """

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
