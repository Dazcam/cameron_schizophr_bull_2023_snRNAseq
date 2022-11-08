# -------------------------------------------------------------------------------------
#
#
#    Script for running LDSC on snRNA-seq data
#
#
# -------------------------------------------------------------------------------------

# ---------  SET SMK PARAMS  ----------
configfile: "../config/config.yaml"

# -------------  RULES  ---------------
rule ldsr_make_annot:
    # Input can be bed file with gene boundaries or gene set with separate gene coord file
    input:   gene_set = "../results/gene_lists/shi_bc/LDSR/{CELL_TYPE}.{GENE_WINDOW}.bed",
             bim_file = "../resources/ldsr/reference_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.{CHR}.bim"
    output:  "../results/LDSR_annotation_files/snRNAseq.{CELL_TYPE}.{GENE_WINDOW}.{CHR}.annot.gz"
    conda:   "../envs/ldsr.yml"
    message: "Creating annotation files for snRNAseq: {wildcards.CELL_TYPE}, {wildcards.GENE_WINDOW}, Chr {wildcards.CHR}"
    log:     "../results/logs/ldsc/make_annot.snRNAseq.{CELL_TYPE}.{GENE_WINDOW}.Chr{CHR}.log"
    shell:
             """
             WINDOW=({wildcards.GENE_WINDOW})
             if [ "${{WINDOW}}" == "0UP_0DOWN" ]; then
             
             echo $WINDOW

             python ../resources/ldsr/make_annot.py \
             --bed-file {input.gene_set} \
             --windowsize 100000 \
             --bimfile {input.bim_file} \
             --annot-file {output} 2> {log}
             
             else 
          
             echo $WINDOW
             
             python ../resources/ldsr/make_annot.py \
             --bed-file {input.gene_set} \
             --windowsize 100000 \
             --bimfile {input.bim_file} \
             --annot-file {output} 2> {log}             


             fi
             """
        
rule ldsr_ld_scores:
    input:   annot = "../results/LDSR_annotation_files/snRNAseq.{CELL_TYPE}.{GENE_WINDOW}.{CHR}.annot.gz",
             bfile_folder = "../resources/ldsr/reference_files/1000G_EUR_Phase3_plink",
             snps_folder = "../resources/ldsr/reference_files/hapmap3_snps"
    output:  "../results/LDSR_annotation_files/snRNAseq.{CELL_TYPE}.{GENE_WINDOW}.{CHR}.l2.ldscore.gz"
    conda:   "../envs/ldsr.yml"
    params:  bfile = "../resources/ldsr/reference_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.{CHR}",
             ldscores = "../results/LDSR_annotation_files/snRNAseq.{CELL_TYPE}.{GENE_WINDOW}.{CHR}",
             snps = "../resources/ldsr/reference_files/w_hm3.snplist_rsIds"
    message: "Running LDSR Phase 3 for {wildcards.CELL_TYPE}, {wildcards.GENE_WINDOW}, CHR {wildcards.CHR}" 
    log:     "../results/logs/LDSR/snRNAseq.{CELL_TYPE}.{GENE_WINDOW}.Chr{CHR}_ldsc.log"
    shell:
        "python ../resources/ldsr/ldsc.py --thin-annot --l2 --bfile {params.bfile} --ld-wind-cm 1 "
        "--annot {input.annot} --out {params.ldscores} --print-snps {params.snps} 2> {log}"


rule ldsr_stratified_baseline_v12:
    input:   GWAS = "../results/GWAS_for_ldsr/{GWAS}_hg19_ldsc_ready.sumstats.gz",
             LDSR = expand("../results/LDSR_annotation_files/snRNAseq.{CELL_TYPE}.{GENE_WINDOW}.{CHR}.l2.ldscore.gz", CELL_TYPE = config["RNA_CELL_TYPES"], GENE_WINDOW = config["GENE_WINDOW"], CHR = range(1,23))
    output:  "../results/LDSR_part_herit/baseline_v1.2/snRNAseq.{CELL_TYPE}.{GENE_WINDOW}.{GWAS}_baseline.v1.2.results"
    conda:   "../envs/ldsr.yml"
    params:  weights = "../resources/ldsr/reference_files/weights_hm3_no_hla/weights.",
             baseline = "../resources/ldsr/reference_files/baseline_v1.2_1000G_Phase3/baseline.",
             frqfile = "../resources/ldsr/reference_files/1000G_Phase3_frq/1000G.EUR.QC.",
             LD_anns = "../results/LDSR_annotation_files/snRNAseq.{CELL_TYPE}.{GENE_WINDOW}.",
             out_file = "../results/LDSR_part_herit/baseline_v1.2/snRNAseq.{CELL_TYPE}.{GENE_WINDOW}.{GWAS}_baseline.v1.2"
    message: "Running Prt Hrt with {wildcards.CELL_TYPE} {wildcards.GENE_WINDOW} and {wildcards.GWAS} GWAS"
    log:     "../results/logs/LDSR/snRNAseq.{CELL_TYPE}.{GENE_WINDOW}.{GWAS}.baseline.v1.2_partHerit.log"
    shell:
             "python ../resources/ldsr/ldsc.py --h2 {input.GWAS} --w-ld-chr {params.weights} "
             "--ref-ld-chr {params.baseline},{params.LD_anns} --overlap-annot "
             "--frqfile-chr {params.frqfile} --out {params.out_file} --print-coefficients 2> {log}"


rule ldsr_stratified_summary:
    # This is still optimised for multiple quantiles so creating > 100 single line files
    input:   expand("../results/LDSR_part_herit/baseline_v1.2/snRNAseq.{CELL_TYPE}.{GENE_WINDOW}.{GWAS}_baseline.v1.2.results", CELL_TYPE = config["RNA_CELL_TYPES"], GENE_WINDOW = config["GENE_WINDOW"], GWAS = config["GWAS"])
    output:  "../results/LDSR_part_herit/baseline_v1.2/snRNAseq_LDSR_{GWAS}_baseline.v1.2_summary.tsv"
    message: "Creating summary file for {wildcards.GWAS} GWAS"
    params:  dir = "../results/LDSR_part_herit/baseline_v1.2/",
             cell_types = "../resources/sheets/GE_celltypes.tsv"
    log:     "../results/logs/LDSR/snRNAseq.{GWAS}_baseline.v1.2_partHerit.summary.log"
    shell:
             """

             
             head -1 {params.dir}snRNAseq.LGE.10UP_10DOWN.SCZ_baseline.v1.2.results > {output}
             File={params.cell_types}
             Lines=$(cat $File)
             for Line in $Lines
             do
             grep L2_1 ../results/LDSR_part_herit/baseline_v1.2/snRNAseq."$Line".{wildcards.GWAS}_baseline.v1.2.results | sed "s/L2_1/$Line/g" >> {output} 2> {log}
             done

             """

#rule create_top_decile_tables:
#    input:   expand("../results/LDSR_part_herit/baseline_v1.2/snRNAseq_LDSC_{CELL_TYPE}_{GWAS}_baseline.v1.2_summary.tsv", CELL_TYPE = config["RNA_CELL_TYPES"], GWAS = config["SUMSTATS"])
#    output:  "../results/LDSR_part_herit/baseline_v1.2/snRNAseq_LDSC_{GWAS}_baseline.v1.2_top10pc.tsv"
#    message: "Creating LDSC top decile tables for {wildcards.GWAS} GWAS"
#    params:  dir = "../results/LDSR_part_herit/baseline_v1.2/"
#    log:     "../results/logs/LDSR/snRNAseq.{GWAS}_partHerit_baseline.v1.2_top10pc_summary.log"
#    shell:
#             """
#             head -1 {params.dir}snRNAseq_LDSC_Cer-OPC_Q10_SCZ_baseline.v1.2.results > {output}
#             for file in `ls {params.dir}*Q10_{wildcards.GWAS}*`; do
#             CELL_TYPE=$(echo ${{file}} | cut -d'_' -f6)
#             tail -1 ${{file}} >> {output}
#             sed -i "s/L2_1/${{CELL_TYPE}}/g" {output}
#             sed -i '/Total time elapsed/d' {output}
#             done
#             """        
        
        
# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
