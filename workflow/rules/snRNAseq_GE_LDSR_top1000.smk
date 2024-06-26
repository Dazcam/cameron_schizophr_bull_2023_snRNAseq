# -------------------------------------------------------------------------------------
#
#
#    Script for running LDSC on snRNA-seq  data top 1000 GeX specifcity genes
#
#
# -------------------------------------------------------------------------------------

# ---------  SET SMK PARAMS  ----------
configfile: "../config/config.yaml"

# -------------  RULES  ---------------
rule ldsr_make_annot_top1000:
    # Input can be bed file with gene boundaries or gene set with separate gene coord file
    input:   gene_set = "../results/02GENE_LISTS/shi_bc/LDSR_top_1000_genes/{CELL_TYPE}.{GENE_WINDOW}.bed",
             bim_file = "../resources/ldsr/reference_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.{CHR}.bim"
    output:  "../results/05LDSR/annotation_files/LDSR_top_1000_genes/snRNAseq.{CELL_TYPE}.{GENE_WINDOW}.{CHR}.annot.gz"
    conda:   "../envs/ldsr.yml"
    message: "Creating annotation files for snRNAseq top 1000 genes/: {wildcards.CELL_TYPE}, {wildcards.GENE_WINDOW}, Chr {wildcards.CHR}"
    log:     "../results/00LOG/05LDSR/LDSR_top_1000_genes/make_annot.snRNAseq.{CELL_TYPE}.{GENE_WINDOW}.Chr{CHR}.log"
    shell:
             """

             python ../resources/ldsr/make_annot.py \
             --bed-file {input.gene_set} \
             --windowsize 0 \
             --bimfile {input.bim_file} \
             --annot-file {output} 2> {log}

             """

rule ldsr_ld_scores_top1000:
    input:   annot = "../results/05LDSR/annotation_files/LDSR_top_1000_genes/snRNAseq.{CELL_TYPE}.{GENE_WINDOW}.{CHR}.annot.gz",
             bfile_folder = "../resources/ldsr/reference_files/1000G_EUR_Phase3_plink",
             snps_folder = "../resources/ldsr/reference_files/hapmap3_snps"
    output:  "../results/05LDSR/annotation_files/LDSR_top_1000_genes/snRNAseq.{CELL_TYPE}.{GENE_WINDOW}.{CHR}.l2.ldscore.gz"
    conda:   "../envs/ldsr.yml"
    params:  bfile = "../resources/ldsr/reference_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.{CHR}",
             ldscores = "../results/05LDSR/annotation_files/LDSR_top_1000_genes/snRNAseq.{CELL_TYPE}.{GENE_WINDOW}.{CHR}",
             snps = "../resources/ldsr/reference_files/w_hm3.snplist_rsIds"
    message: "Running LDSR Phase 3 for {wildcards.CELL_TYPE}, top 1000 genes {wildcards.GENE_WINDOW}, CHR {wildcards.CHR}"
    log:     "../results/00LOG/05LDSR/LDSR_top_1000_genes/snRNAseq.{CELL_TYPE}.{GENE_WINDOW}.Chr{CHR}_ldsc.log"
    shell:
        "python ../resources/ldsr/ldsc.py --thin-annot --l2 --bfile {params.bfile} --ld-wind-cm 1 "
        "--annot {input.annot} --out {params.ldscores} --print-snps {params.snps} 2> {log}"


rule ldsr_stratified_baseline_v12_top1000:
    input:   GWAS = "../results/03SUMSTATS/{GWAS}_hg19_LDSR_ready.sumstats.gz",
             LDSR = expand("../results/05LDSR/annotation_files/LDSR_top_1000_genes/snRNAseq.{CELL_TYPE}.{GENE_WINDOW}.{CHR}.l2.ldscore.gz", CELL_TYPE = config["CELL_TYPES"], GENE_WINDOW = config["GENE_WINDOW"], CHR = range(1,23))
    output:  "../results/05LDSR/part_herit/baseline_v1.2/LDSR_top_1000_genes/snRNAseq.{CELL_TYPE}.{GENE_WINDOW}.{GWAS}_baseline.v1.2.results"
    conda:   "../envs/ldsr.yml"
    params:  weights = "../resources/ldsr/reference_files/weights_hm3_no_hla/weights.",
             baseline = "../resources/ldsr/reference_files/baseline_v1.2_1000G_Phase3/baseline.",
             frqfile = "../resources/ldsr/reference_files/1000G_Phase3_frq/1000G.EUR.QC.",
             LD_anns = "../results/05LDSR/annotation_files/LDSR_top_1000_genes/snRNAseq.{CELL_TYPE}.{GENE_WINDOW}.",
             out_file = "../results/05LDSR/part_herit/baseline_v1.2/LDSR_top_1000_genes/snRNAseq.{CELL_TYPE}.{GENE_WINDOW}.{GWAS}_baseline.v1.2"
    message: "Running Prt Hrt with {wildcards.CELL_TYPE} top 1000 genes {wildcards.GENE_WINDOW} and {wildcards.GWAS} GWAS"
    log:     "../results/00LOG/05LDSR/LDSR_top_1000_genes/snRNAseq.{CELL_TYPE}.{GENE_WINDOW}.{GWAS}.baseline.v1.2_partHerit.log"
    shell:
             "python ../resources/ldsr/ldsc.py --h2 {input.GWAS} --w-ld-chr {params.weights} "
             "--ref-ld-chr {params.baseline},{params.LD_anns} --overlap-annot "
             "--frqfile-chr {params.frqfile} --out {params.out_file} --print-coefficients 2> {log}"


rule ldsr_stratified_summary_top1000:
    input:   expand("../results/05LDSR/part_herit/baseline_v1.2/LDSR_top_1000_genes/snRNAseq.{CELL_TYPE}.{GENE_WINDOW}.{GWAS}_baseline.v1.2.results", CELL_TYPE = config["CELL_TYPES"], GENE_WINDOW = config["GENE_WINDOW"], GWAS = config["GWAS"])
    output:  "../results/05LDSR/part_herit/baseline_v1.2/LDSR_top_1000_genes/snRNAseq_LDSR_{GWAS}_baseline.v1.2_summary.tsv"
    message: "Creating summary file for {wildcards.GWAS} GWAS top 1000 genes"
    params:  dir = "../results/05LDSR/part_herit/baseline_v1.2/LDSR_top_1000_genes/",
             cell_types = "../resources/sheets/GE_celltypes.tsv"
    log:     "../results/00LOG/05LDSR/LDSR_top_1000_genes/snRNAseq.{GWAS}_baseline.v1.2_partHerit.summary.log"
    shell:
             """


             head -1 {params.dir}snRNAseq.LGE.100UP_100DOWN.SCZ_baseline.v1.2.results > {output}
             File={params.cell_types}
             Lines=$(cat $File)
             for Line in $Lines
             do
             grep L2_1 ../results/05LDSR/part_herit/baseline_v1.2/LDSR_top_1000_genes/snRNAseq."$Line".{wildcards.GWAS}_baseline.v1.2.results | sed "s/L2_1/$Line/g" >> {output} 2> {log}
             done

             """


# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
