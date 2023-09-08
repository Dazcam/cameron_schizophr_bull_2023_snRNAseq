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
    log:     "../results/01LOG/03MAGMA/snRNAseq.GE.getRefs.log"
    shell:
             """
             
             wget https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip -P {params}
             unzip {params}g1000_eur.zip             

             """        

rule magma_map_snps_to_genes:
    input:   snp_loc = "../resources/refs/g1000_eur.bim",
             gene_loc = "../resources/refs/NCBI37.3.MHCremoved.gene.loc.txt"
    output:  "../results/03MAGMA/snRNAseq_GE.magma.{GENE_WINDOW}.genes.annot"
    params:  "../results/03MAGMA/snRNAseq_GE.magma.{GENE_WINDOW}"
    message: "Running MAGMA annotation step to map SNPs to genes. Gene window: {wildcards.GENE_WINDOW}"
    log:     "../results/00LOG/03MAGMA/snRNAseq.GE.annotate.snps2genes.{GENE_WINDOW}.log"
    run: 
        if "0UP_0DOWN" in wildcards.GENE_WINDOW:

            print("\nMap SNPs to genes for gene window: ", wildcards.GENE_WINDOW, "\n")

            shell("""

            module load magma/1.10
            magma --annotate window=0,0 --snp-loc {input.snp_loc} --gene-loc {input.gene_loc} --out {params} &> {log}

            """)

        elif "10UP_10DOWN" in wildcards.GENE_WINDOW:
        
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
    input:   gene_annot = "../results/03MAGMA/snRNAseq_GE.magma.{GENE_WINDOW}.genes.annot",
             gwas = "../results/GWAS_for_MAGMA/{GWAS}_hg19_magma_ready.tsv"
    output:  "../results/03MAGMA/snRNAseq_GE_{GWAS}.magma.{GENE_WINDOW}.genes.raw",
             "../results/03MAGMA/snRNAseq_GE_{GWAS}.magma.{GENE_WINDOW}.genes.out"
    params:  ref = "../resources/refs/g1000_eur",
             out = "../results/03MAGMA/snRNAseq_GE_{GWAS}.magma.{GENE_WINDOW}"
    message: "Running magma gene analysis step for {wildcards.GWAS}. Gene window: {wildcards.GENE_WINDOW}"
    log:     "../results/00LOG/03MAGMA/snRNAseq.GE.gene_analysis.{GWAS}.{GENE_WINDOW}.log"
    shell:
             """

             module load magma/1.10
             magma --bfile {params.ref} --pval {input.gwas} ncol='N' --gene-annot {input.gene_annot} --out {params.out} &> {log}

             """

rule magma_gene_set_analysis:
    input:   genes = "../results/03MAGMA/snRNAseq_GE_{GWAS}.magma.{GENE_WINDOW}.genes.raw",
             data  = "../results/02GENE_LISTS/{SEURAT_OBJ}/MAGMA/shi_top10_lvl_{LEVEL}.txt"
    output:  "../results/03MAGMA/snRNAseq_GE_{GWAS}.{SEURAT_OBJ}.lvl_{LEVEL}.magma.{GENE_WINDOW}.gsa.out"
    params:  out = "../results/03MAGMA/snRNAseq_GE_{GWAS}.{SEURAT_OBJ}.lvl_{LEVEL}.magma.{GENE_WINDOW}"
    message: "Running MAGMA gene set analysis step for {wildcards.GWAS}, {wildcards.SEURAT_OBJ}, cluster level {wildcards.LEVEL}. Gene window: {wildcards.GENE_WINDOW}"
    log:     "../results/00LOG/03MAGMA/snRNAseq.GE.gene_set_analysis.{GWAS}.{SEURAT_OBJ}.lvl_{LEVEL}.{GENE_WINDOW}.log"
    shell:
             """

             module load magma/1.10
             magma --gene-results {input.genes} --set-annot {input.data} --out {params.out} &> {log}

             """

rule magma_conditional:
    input:   gene_list = "../results/02GENE_LISTS/shi_bc/MAGMA_CONDITIONAL/skene_bryois_InN_entrez_gene_list.tsv",
             scz_magma = "../results/03MAGMA/snRNAseq_GE_SCZ.magma.{GENE_WINDOW}.genes.raw" 
    output:  "../results/03MAGMA/magma_conditional/magma_all_sig_and_skene_condition_{CONDITION}.{GENE_WINDOW}.gsa.out"
    params:  "../results/03MAGMA/magma_conditional/magma_all_sig_and_skene_condition_{CONDITION}.{GENE_WINDOW}"
    message: "Running MAGMA on all significant cell types conditioning on {wildcards.CONDITION}. Gene window: {wildcards.GENE_WINDOW}"
    log:     "../results/00LOG/03MAGMA/magma_conditional/snRNAseq.magma.conditional.{CONDITION}.{GENE_WINDOW}.log"
    shell:
             """

             magma --gene-results {input.scz_magma} --set-annot {input.gene_list} --model condition={wildcards.CONDITION} --out {params} &> {log}

             """

rule magma_gene_set_analysis_top1000:
    input:   genes = "../results/03MAGMA/snRNAseq_GE_{GWAS}.magma.{GENE_WINDOW}.genes.raw",
             data  = "../results/02GENE_LISTS/{SEURAT_OBJ}/MAGMA/shi_top1000_lvl_{LEVEL}.txt"
    output:  "../results/03MAGMA/snRNAseq_GE_{GWAS}.{SEURAT_OBJ}.top1000.lvl_{LEVEL}.magma.{GENE_WINDOW}.gsa.out"
    params:  out = "../results/03MAGMA/snRNAseq_GE_{GWAS}.{SEURAT_OBJ}.top1000.lvl_{LEVEL}.magma.{GENE_WINDOW}"
    message: "Running MAGMA gene set analysis step for {wildcards.GWAS}, {wildcards.SEURAT_OBJ}, top 1000 genes, cluster level {wildcards.LEVEL}. Gene window: {wildcards.GENE_WINDOW}"
    log:     "../results/00LOG/03MAGMA/snRNAseq.GE.gene_set_analysis.{GWAS}.{SEURAT_OBJ}.top1000.lvl_{LEVEL}.{GENE_WINDOW}.log"
    shell:
             """

             module load magma/1.10
             magma --gene-results {input.genes} --set-annot {input.data} --out {params.out} &> {log}

             """

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
