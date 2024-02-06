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
    log:     "../results/01LOG/04MAGMA/snRNAseq.GE.getRefs.log"
    shell:
             """
             
             wget https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip -P {params}
             unzip {params}g1000_eur.zip             

             """        

rule magma_map_snps_to_genes:
    input:   snp_loc = "../resources/refs/g1000_eur.bim",
             gene_loc = "../resources/refs/NCBI37.3.MHCremoved.gene.loc.txt"
    output:  "../results/04MAGMA/snRNAseq_GE.magma.{GENE_WINDOW}.genes.annot"
    params:  "../results/04MAGMA/snRNAseq_GE.magma.{GENE_WINDOW}"
    message: "Running MAGMA annotation step to map SNPs to genes. Gene window: {wildcards.GENE_WINDOW}"
    log:     "../results/00LOG/04MAGMA/snRNAseq.GE.annotate.snps2genes.{GENE_WINDOW}.log"
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
    input:   gene_annot = "../results/04MAGMA/snRNAseq_GE.magma.{GENE_WINDOW}.genes.annot",
             gwas = "../results/03SUMSTATS/{GWAS}_hg19_MAGMA_ready.tsv"
    output:  "../results/04MAGMA/snRNAseq_GE_{GWAS}.magma.{GENE_WINDOW}.genes.raw",
             "../results/04MAGMA/snRNAseq_GE_{GWAS}.magma.{GENE_WINDOW}.genes.out"
    params:  ref = "../resources/refs/g1000_eur",
             out = "../results/04MAGMA/snRNAseq_GE_{GWAS}.magma.{GENE_WINDOW}"
    message: "Running magma gene analysis step for {wildcards.GWAS}. Gene window: {wildcards.GENE_WINDOW}"
    log:     "../results/00LOG/04MAGMA/snRNAseq.GE.gene_analysis.{GWAS}.{GENE_WINDOW}.log"
    shell:
             """

             module load magma/1.10
             magma --bfile {params.ref} --pval {input.gwas} ncol='N' --gene-annot {input.gene_annot} --out {params.out} &> {log}

             """

rule magma_gene_set_analysis:
    input:   genes = "../results/04MAGMA/snRNAseq_GE_{GWAS}.magma.{GENE_WINDOW}.genes.raw",
             data  = "../results/02GENE_LISTS/{SEURAT_OBJ}/MAGMA/shi_top10_lvl_{LEVEL}.txt"
    output:  "../results/04MAGMA/snRNAseq_GE_{GWAS}.{SEURAT_OBJ}.lvl_{LEVEL}.magma.{GENE_WINDOW}.gsa.out"
    params:  out = "../results/04MAGMA/snRNAseq_GE_{GWAS}.{SEURAT_OBJ}.lvl_{LEVEL}.magma.{GENE_WINDOW}"
    message: "Running MAGMA gene set analysis step for {wildcards.GWAS}, {wildcards.SEURAT_OBJ}, cluster level {wildcards.LEVEL}. Gene window: {wildcards.GENE_WINDOW}"
    log:     "../results/00LOG/04MAGMA/snRNAseq.GE.gene_set_analysis.{GWAS}.{SEURAT_OBJ}.lvl_{LEVEL}.{GENE_WINDOW}.log"
    shell:
             """

             module load magma/1.10
             magma --gene-results {input.genes} --set-annot {input.data} --out {params.out} &> {log}

             """

rule magma_conditional:
    input:   gene_list = "../results/02GENE_LISTS/shi_bc/MAGMA_CONDITIONAL/snRNAseq_GE_conditional_gene_sets.txt",
             scz_magma = "../results/04MAGMA/snRNAseq_GE_SCZ.magma.{GENE_WINDOW}.genes.raw" 
    output:  "../results/04MAGMA/magma_conditional/snRNAseq_GE_SCZ_magma_conditional_{CONDITION}.{GENE_WINDOW}.gsa.out"
    params:  "../results/04MAGMA/magma_conditional/snRNAseq_GE_SCZ_magma_conditional_{CONDITION}.{GENE_WINDOW}"
    message: "Running MAGMA conditional analyses: Conditioning on {wildcards.CONDITION}. Gene window: {wildcards.GENE_WINDOW}"
    log:     "../results/00LOG/04MAGMA/magma_conditional/snRNAseq.magma.conditional.{CONDITION}.{GENE_WINDOW}.log"
    shell:
             """

             magma --gene-results {input.scz_magma} --set-annot {input.gene_list} --model condition={wildcards.CONDITION} --out {params} &> {log}

             """

rule magma_gene_set_analysis_top1000:
    input:   genes = "../results/03MAGMA/snRNAseq_GE_{GWAS}.magma.{GENE_WINDOW}.genes.raw",
             data  = "../results/02GENE_LISTS/{SEURAT_OBJ}/MAGMA/shi_top1000_lvl_{LEVEL}.txt"
    output:  "../results/03MAGMA/snRNAseq_GE_{GWAS}.{SEURAT_OBJ}.top1000.lvl_{LEVEL}.magma.{GENE_WINDOW}.gsa.out"
    params:  "../results/03MAGMA/snRNAseq_GE_{GWAS}.{SEURAT_OBJ}.top1000.lvl_{LEVEL}.magma.{GENE_WINDOW}"
    message: "Running MAGMA gene set analysis step for {wildcards.GWAS}, {wildcards.SEURAT_OBJ}, top 1000 genes, cluster level {wildcards.LEVEL}. Gene window: {wildcards.GENE_WINDOW}"
    log:     "../results/00LOG/03MAGMA/snRNAseq.GE.gene_set_analysis.{GWAS}.{SEURAT_OBJ}.top1000.lvl_{LEVEL}.{GENE_WINDOW}.log"
    shell:
             """

             module load magma/1.10
             magma --gene-results {input.genes} --set-annot {input.data} --out {params} &> {log}

             """


rule magma_gene_set_analysis_diffExp:
    input:   genes = "../results/03MAGMA/snRNAseq_GE_SCZ.magma.35UP_10DOWN.genes.raw",
             data  = "../results/02GENE_LISTS/shi_bc/MAGMA_DIFFEXP/snRNAseq_GE_diffexp_gene_sets.txt"
    output:  "../results/03MAGMA/magma_diffexp/snRNAseq_GE_magma_diffexp.35UP_10DOWN.gsa.out"
    params:  "../results/03MAGMA/magma_diffexp/snRNAseq_GE_magma_diffexp.35UP_10DOWN"
    message: "Running MAGMA gene set analysis step for SCZ, monocle diff exp genes, Gene window: 35UP_10DOWN"
    log:     "../results/00LOG/03MAGMA/snRNAseq.GE.gene_set_analysis.SCZ.diffexp.35UP_10DOWN.log"
    shell:
             """

             module load magma/1.10
             magma --gene-results {input.genes} --set-annot {input.data} --out {params} &> {log}

             """

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
