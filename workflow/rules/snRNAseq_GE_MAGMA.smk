# -------------------------------------------------------------------------------------
#
#
#    MAGMA analyses
#
#
# -------------------------------------------------------------------------------------


# ---------  SET SMK PARAMS  ----------
configfile: "../config/config.yaml"
localrules: magma_download_refs

# -------------  RULES  ---------------

rule magma_download_refs:
    output:  "../resources/refs/g1000_eur.bim"
    params:  "../resources/refs/"
    log:     "../results/logs/magma/snRNAseq.GE.getRefs.log"
    shell:
             """
             
             wget https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip -P {params}
             unzip {params}g1000_eur.zip             

             """        

rule magma_map_snps_to_genes:
    input:   snp_loc = "../resources/refs/g1000_eur.bim",
             gene_loc = "../resources/refs/NCBI37.3.gene.loc.extendedMHCexcluded.txt"
    output:  "../results/magma/snRNAseq_GE.magma.{GENE_WINDOW}.genes.annot"
    params:  "../results/magma/snRNAseq_GE.magma.{GENE_WINDOW}"
    message: "Running magma annotation step to map SNPs to genes. Gene window: {wildcards.GENE_WINDOW}"
    log:     "../results/logs/magma/snRNAseq.GE.annotate.snps2genes.{GENE_WINDOW}.log"
    run: 
        if "10UP_10DOWN" in wildcards.GENE_WINDOW:
        
            print("\nMap SNPs to genes for gene window: ", wildcards.GENE_WINDOW, "\n")
                           
            shell("""
                 
            module load magma/1.10
            magma --annotate window=10,10 --snp-loc {input.snp_loc} --gene-loc {input.gene_loc} --out {params} &> {log}
                 
            """)
        
        elif "35UP_10DOWN" in wildcards.GENE_WINDOW:
        
            print("\nMap SNPs to genes for gene window: ", wildcards.GENE_WINDOW, "\n")
                           
            shell("""
                 
            module load magma/1.10
            magma --annotate window=35,10 --snp-loc {input.snp_loc} --gene-loc {input.gene_loc} --out {params} &> {log}
                 
            """)

        else:

            print("\nMap SNPs to genes for gene window: ", wildcards.GENE_WINDOW, "\n")

            shell("""

            module load magma/1.10
            magma --annotate window=100,100 --snp-loc {input.snp_loc} --gene-loc {input.gene_loc} --out {params} &> {log}

            """)

rule magma_gene_analysis:
    input:   gene_annot = "../results/magma/snRNAseq_GE.magma.{GENE_WINDOW}.genes.annot",
             gwas = "../results/GWAS_for_MAGMA/{GWAS}_hg19_magma_ready.tsv"
    output:  "../results/magma/snRNAseq_GE_{GWAS}.magma.{GENE_WINDOW}.genes.raw",
             "../results/magma/snRNAseq_GE_{GWAS}.magma.{GENE_WINDOW}.genes.out"
    params:  ref = "../resources/refs/g1000_eur",
             out = "../results/magma/snRNAseq_GE_{GWAS}.magma.{GENE_WINDOW}"
    message: "Running magma gene analysis step for {wildcards.GWAS}. Gene window: {wildcards.GENE_WINDOW}"
    log:     "../results/logs/magma/snRNAseq.GE.gene_analysis.{GWAS}.{GENE_WINDOW}.log"
    shell:
             """

             module load magma/1.10
             magma --bfile {params.ref} --pval {input.gwas} ncol='N' --gene-annot {input.gene_annot} --out {params.out} &> {log}

             """

rule magma_gene_set_analysis:
    input:   genes = "../results/magma/snRNAseq_GE_{GWAS}.magma.{GENE_WINDOW}.genes.raw",
             data  = "../results/gene_lists/{SEURAT_OBJ}/MAGMA/shi_top10_lvl_{LEVEL}.txt"
    output:  "../results/magma/snRNAseq_GE_{GWAS}.{SEURAT_OBJ}.lvl_{LEVEL}.magma.{GENE_WINDOW}.gsa.out"
    params:  out = "../results/magma/snRNAseq_GE_{GWAS}.{SEURAT_OBJ}.lvl_{LEVEL}.magma.{GENE_WINDOW}"
    message: "Running magma gene set analysis step for {wildcards.GWAS}, {wildcards.SEURAT_OBJ}, cluster level {wildcards.LEVEL}. Gene window: {wildcards.GENE_WINDOW}"
    log:     "../results/logs/magma/snRNAseq.GE.gene_set_analysis.{GWAS}.{SEURAT_OBJ}.lvl_{LEVEL}.{GENE_WINDOW}.log"
    shell:
             """

             module load magma/1.10
             magma --gene-results {input.genes} --set-annot {input.data} --out {params.out} &> {log}

             """

rule magma_conditional:
    input:   gene_list = "../results/gene_lists/MAGMA/shi_top10_lvl_2_sig_plus_skene.txt",
             scz_magma = "../results/magma/snRNAseq_GE_SCZ.magma.{GENE_WINDOW}.genes.raw" 
    output:  "../results/magma_conditional/magma_all_sig_and_skene_condition_{CONDITION}.{GENE_WINDOW}.gsa.out"
    params:  "../results/magma_conditional/magma_all_sig_and_skene_condition_{CONDITION}.{GENE_WINDOW}"
    message: "Running magma on all significant cell types conditioning on {wildcards.CONDITION}. Gene window: {wildcards.GENE_WINDOW}"
    log:     "../results/logs/magma_conditional/snRNAseq.magma.conditional.{CONDITION}.{GENE_WINDOW}.log"
    shell:
             """

             magma --gene-results {input.scz_magma} --set-annot {input.gene_list} --model condition={wildcards.CONDITION} --out {params} &> {log}

             """

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
