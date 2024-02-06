# -------------------------------------------------------------------------------------
#
#
#    Script for running LDSC on snRNA-seq downsampled data
#
#
# -------------------------------------------------------------------------------------

# ---------  SET SMK PARAMS  ----------
configfile: "../config/config.yaml"

# -------------  RULES  ---------------
rule ldsr_make_annot_cond_public:
    # Input can be bed file with gene boundaries or gene set with separate gene coord file
    input:   gene_set = "../results/02GENE_LISTS/shi_bc/LDSR/{PUBLIC_CELL_TYPE}.100UP_100DOWN.bed",
             bim_file = "../resources/ldsr/reference_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.{CHR}.bim"
    output:  "../results/05LDSR/annotation_files/LDSR_cond_public/snRNAseq.{PUBLIC_CELL_TYPE}.100UP_100DOWN.{CHR}.annot.gz"
    conda:   "../envs/ldsr.yml"
    message: "Creating annotation files for snRNAseq public conditional analysis: {wildcards.PUBLIC_CELL_TYPE}, 100UP_100DOWN, Chr {wildcards.CHR}"
    log:     "../results/00LOG/05LDSR/shi_bc/make_annot.snRNAseq.{PUBLIC_CELL_TYPE}.100UP_100DOWN.Chr{CHR}.log"
    shell:
             """

             python ../resources/ldsr/make_annot.py \
             --bed-file {input.gene_set} \
             --windowsize 0 \
             --bimfile {input.bim_file} \
             --annot-file {output} 2> {log}

             """

rule ldsr_ld_scores_cond_public:
    input:   annot = "../results/05LDSR/annotation_files/LDSR_cond_public/snRNAseq.{PUBLIC_CELL_TYPE}.100UP_100DOWN.{CHR}.annot.gz",
             bfile_folder = "../resources/ldsr/reference_files/1000G_EUR_Phase3_plink",
             snps_folder = "../resources/ldsr/reference_files/hapmap3_snps"
    output:  "../results/05LDSR/annotation_files/LDSR_cond_public/snRNAseq.{PUBLIC_CELL_TYPE}.100UP_100DOWN.{CHR}.l2.ldscore.gz"
    conda:   "../envs/ldsr.yml"
    params:  bfile = "../resources/ldsr/reference_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.{CHR}",
             ldscores = "../results/05LDSR/annotation_files/LDSR_cond_public/snRNAseq.{PUBLIC_CELL_TYPE}.100UP_100DOWN.{CHR}",
             snps = "../resources/ldsr/reference_files/w_hm3.snplist_rsIds"
    message: "Running LDSR Phase 3 for {wildcards.PUBLIC_CELL_TYPE}, 100UP_100DOWN, CHR {wildcards.CHR}"
    log:     "../results/00LOG/05LDSR/LDSR_cond_public/snRNAseq.{PUBLIC_CELL_TYPE}.100UP_100DOWN.Chr{CHR}_ldsc.log"
    shell:
	     "python ../resources/ldsr/ldsc.py --thin-annot --l2 --bfile {params.bfile} --ld-wind-cm 1 "
             "--annot {input.annot} --out {params.ldscores} --print-snps {params.snps} 2> {log}"

rule ldsr_stratified_baseline_v12_cond_public:
    input:   GWAS = "../results/03SUMSTATS/{GWAS}_hg19_LDSR_ready.sumstats.gz",
             LDSR = expand("../results/05LDSR/annotation_files/snRNAseq.{COND_CELL_A}.100UP_100DOWN.{CHR}.l2.ldscore.gz", COND_CELL_A = config["COND_CELL_A"], CHR = range(1,23)),
             COND = expand("../results/05LDSR/annotation_files/LDSR_cond_public/snRNAseq.{PUBLIC_CELL_TYPE}.100UP_100DOWN.{CHR}.l2.ldscore.gz", PUBLIC_CELL_TYPE = config["PUBLIC_CELL_TYPES"], CHR = range(1,23))
    output:  "../results/05LDSR/part_herit/baseline_v1.2/LDSR_cond_public/snRNAseq.{COND_CELL_A}_vs_{PUBLIC_CELL_TYPE}.100UP_100DOWN.{GWAS}_baseline.v1.2.results"
    conda:   "../envs/ldsr.yml"
    params:  weights = "../resources/ldsr/reference_files/weights_hm3_no_hla/weights.",
             baseline = "../resources/ldsr/reference_files/baseline_v1.2_1000G_Phase3/baseline.",
             frqfile = "../resources/ldsr/reference_files/1000G_Phase3_frq/1000G.EUR.QC.",
             LD_anns = "../results/05LDSR/annotation_files/snRNAseq.{COND_CELL_A}.100UP_100DOWN.",
             cond_anns = "../results/05LDSR/annotation_files/LDSR_cond_public/snRNAseq.{PUBLIC_CELL_TYPE}.100UP_100DOWN.",
             out_file = "../results/05LDSR/part_herit/baseline_v1.2/LDSR_cond_public/snRNAseq.{COND_CELL_A}_vs_{PUBLIC_CELL_TYPE}.100UP_100DOWN.{GWAS}_baseline.v1.2"
    message: "Running Prt Hrt with {wildcards.COND_CELL_A} vs. {wildcards.PUBLIC_CELL_TYPE}, 100UP_100DOWN and {wildcards.GWAS} GWAS"
    log:     "../results/00LOG/05LDSR/LDSR_cond_public/snRNAseq.{COND_CELL_A}_vs_{PUBLIC_CELL_TYPE}.100UP_100DOWN.{GWAS}.baseline.v1.2_partHerit.log"
    shell:
             "python ../resources/ldsr/ldsc.py --h2 {input.GWAS} --w-ld-chr {params.weights} "
             "--ref-ld-chr {params.baseline},{params.cond_anns},{params.LD_anns} --overlap-annot " 
             "--frqfile-chr {params.frqfile} --out {params.out_file} --print-coefficients 2> {log}"

rule ldsr_stratified_summary_cond_public:
    input:   expand("../results/05LDSR/part_herit/baseline_v1.2/LDSR_cond_public/snRNAseq.{COND_CELL_A}_vs_{PUBLIC_CELL_TYPE}.100UP_100DOWN.{GWAS}_baseline.v1.2.results", COND_CELL_A = config['COND_CELL_A'], PUBLIC_CELL_TYPE = config["PUBLIC_CELL_TYPES"], GWAS = config['GWAS'])
    output:  "../results/05LDSR/part_herit/baseline_v1.2/LDSR_cond_public/snRNAseq_LDSR_{GWAS}_baseline.v1.2_summary.tsv"
    message: "Creating adult conditional summary file for {wildcards.GWAS} GWAS"
    params:  dir = "../results/05LDSR/part_herit/baseline_v1.2/LDSR_cond_public/",
             cell_types = "../resources/sheets/GE_celltypes_conditional_public.tsv"
    log:     "../results/00LOG/05LDSR/LDSR_cond_public/snRNAseq.{GWAS}_baseline.v1.2_partHerit.summary.log"
    shell:
             """


             head -1 {params.dir}snRNAseq.CGE_1_vs_mouse_InN.100UP_100DOWN.SCZ_baseline.v1.2.results > {output}
             File={params.cell_types}
             Lines=$(cat $File)
             for Line in $Lines
             do
             grep L2_2 ../results/05LDSR/part_herit/baseline_v1.2/LDSR_cond_public/snRNAseq."$Line".{wildcards.GWAS}_baseline.v1.2.results | sed "s/L2_2/$Line/g" >> {output} 2> {log}
             done

             """


# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
