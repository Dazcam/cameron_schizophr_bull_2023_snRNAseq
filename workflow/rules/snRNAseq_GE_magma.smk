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

rule magma_download_refs:
    output:  dir("../resources/refs/g1000_eur")
    params:  "../resources/refs/"
    log:     "../results/logs/magma/snRNAseq.GE.getRefs.{GWAS}.log"
    shell:
             """
             
             wget https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip -P {params}
             unzip {params}g1000_eur.zip             

             """        

rule magma_map_snps_to_genes:
    input:   snp_loc = 
             gene_loc = "../resources/refs/NCBI37.3.gene.loc.extendedMHCexcluded.txt"
    output:  "../results/magma/snRNAseq_GE_{GWAS}.genes.annot"
    params:  "../results/magma/snRNAseq_GE_{GWAS}"
    message: "Running magma annotation step for {wildcards.GWAS}"
    log:     "../results/logs/magma/snRNAseq.GE.snps2genes.{GWAS}.log"
    shell:
             """

             magma --annotate --window=35,10 --snp-loc {input.snp_loc} --gene-loc {input.gene_loc} --out {params}              

             """

rule magma_gene_analysis:
    input:   gene_annot = "../results/magma/snRNAseq_GE_{GWAS}.genes.annot",
             ref = rules.magma_download_refs.output,
             pval = 
    output:  "../results/magma/snRNAseq_GE_{GWAS}.genes.raw"
    params:  ref = "../resources/refs/g1000_eur"
             out = "../results/magma/snRNAseq_GE_{GWAS}"
    message: "Running magma gene analysis step for {wildcards.GWAS}"
    log:     "../results/logs/magma/snRNAseq.GE.gene_analysis.{GWAS}.log"
    shell:
             """

             magma --bfile {params.ref} --pval [PVAL_FILE] N=[N] --gene-annot {input.gene_annot} --out {params.out}

             """

rule magma_gene_set_analysis:
    input:   genes = "../results/magma/snRNAseq_GE_{GWAS}.genes.raw"
             data  = "../results/gene_lists/MAGMA/shi_top10.txt"
    output:  "../results/magma/snRNAseq_GE_{GWAS}.gsa.out"
    params:  out = "../results/magma/snRNAseq_GE_{GWAS}"
    message: "Running magma gene set analysis step for {wildcards.GWAS}"
    log:     "../results/logs/magma/snRNAseq.GE.gene_analysis.{GWAS}.log"
    shell:
             """

             magma --gene-results {input.genes} --gene-covar {input.data} --out {params.out}

             """

#rule magma_conditional:
#    input:   gene_list = "../results/gene_lists/q10_gene_lists/ALL_SIG_AND_SKENE_entrez_gene_list.tsv",
#             scz_magma = MAGMA_DIR + "SCZ_hg19_magma_ready.sumstats.tsv.10UP.1.5DOWN.genes.raw"
#    output:  "../results/magma_conditional/magma_all_sig_and_skene_condition_{CONDITION}.gsa.out"
#    params:  "../results/magma_conditional/magma_all_sig_and_skene_condition_{CONDITION}"
#    message: "Running magma on all significant cell types conditioning on {wildcards.CONDITION}"
#    log:     "../results/logs/magma_conditional/snRNAseq.magma.conditional.{CONDITION}.log"
#    shell:
#             """

#             magma --gene-results {input.scz_magma} --set-annot {input.gene_list} --model condition={wildcards.CONDITION} --out {params}

#             """

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
