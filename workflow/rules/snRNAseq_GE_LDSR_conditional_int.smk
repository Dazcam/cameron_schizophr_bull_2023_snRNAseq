rule ldsr_stratified_baseline_v12_internal_cond:
    input:   GWAS = "../results/GWAS_for_LDSR/{GWAS}_hg19_ldsc_ready.sumstats.gz",
             LDSR = expand("../results/04LDSR/annotation_files/snRNAseq.{COND_CELL_A}.100UP_100DOWN.{CHR}.l2.ldscore.gz", COND_CELL_A = config["COND_CELL_A"], CHR = range(1,23)),
             COND = expand("../results/04LDSR/annotation_files/snRNAseq.{COND_CELL_B}.100UP_100DOWN.{CHR}.l2.ldscore.gz", COND_CELL_B = config["COND_CELL_B"], CHR = range(1,23))
    output:  "../results/04LDSR/part_herit/baseline_v1.2/LDSR_cond_int/snRNAseq.{COND_CELL_A}_vs_{COND_CELL_B}.100UP_100DOWN.{GWAS}_baseline.v1.2.results"
    conda:   "../envs/ldsr.yml"
    params:  weights = "../resources/ldsr/reference_files/weights_hm3_no_hla/weights.",
             baseline = "../resources/ldsr/reference_files/baseline_v1.2_1000G_Phase3/baseline.",
             frqfile = "../resources/ldsr/reference_files/1000G_Phase3_frq/1000G.EUR.QC.",
             LD_anns = "../results/04LDSR/annotation_files/snRNAseq.{COND_CELL_A}.100UP_100DOWN.",
             cond_anns = "../results/04LDSR/annotation_files/snRNAseq.{COND_CELL_B}.100UP_100DOWN.",
             out_file = "../results/04LDSR/part_herit/baseline_v1.2/LDSR_cond_int/snRNAseq.{COND_CELL_A}_vs_{COND_CELL_B}.100UP_100DOWN.{GWAS}_baseline.v1.2"
    message: "Running Prt Hrt with {wildcards.COND_CELL_A} vs. {wildcards.COND_CELL_B}, 100UP_100DOWN and {wildcards.GWAS} GWAS"
    log:     "../results/00LOG/04LDSR/cond_int/snRNAseq.{COND_CELL_A}_vs_{COND_CELL_B}.100UP_100DOWN.{GWAS}.baseline.v1.2_partHerit.log"
    shell:
             """
             if [[ {wildcards.COND_CELL_A} == {wildcards.COND_CELL_B} ]] 
             then 
               printf 'L2_2\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\n' > {output} 2> {log} 
             else  
               python ../resources/ldsr/ldsc.py --h2 {input.GWAS} --w-ld-chr {params.weights} --ref-ld-chr {params.baseline},{params.cond_anns},{params.LD_anns} --overlap-annot --frqfile-chr {params.frqfile} --out {params.out_file} --print-coefficients 2> {log}
             fi
             """ 

rule ldsr_stratified_summary_internal_cond:
    input:   expand("../results/04LDSR/part_herit/baseline_v1.2/LDSR_cond_int/snRNAseq.{COND_CELL_A}_vs_{COND_CELL_B}.100UP_100DOWN.{GWAS}_baseline.v1.2.results", COND_CELL_A = config['COND_CELL_A'], COND_CELL_B = config['COND_CELL_B'], GWAS = config['GWAS'])
    output:  "../results/04LDSR/part_herit/baseline_v1.2/LDSR_cond_int/snRNAseq_LDSR_{GWAS}_baseline.v1.2_summary.tsv"
    message: "Creating internal conditional summary file for {wildcards.GWAS} GWAS"
    params:  dir = "../results/04LDSR/part_herit/baseline_v1.2/LDSR_cond_int/",
             cell_types = "../resources/sheets/GE_celltypes_conditional_int.tsv"
    log:     "../results/00LOG/04LDSR/cond_int/snRNAseq.{GWAS}_baseline.v1.2_partHerit.summary.log"
    shell:
             """


             head -1 {params.dir}snRNAseq.CGE_1_vs_CGE_2.100UP_100DOWN.SCZ_baseline.v1.2.results > {output}
             File={params.cell_types}
             Lines=$(cat $File)
             for Line in $Lines
             do
             grep L2_2 ../results/04LDSR/part_herit/baseline_v1.2/LDSR_cond_int/snRNAseq."$Line".{wildcards.GWAS}_baseline.v1.2.results | sed "s/L2_2/$Line/g" >> {output} 2> {log}
             done

             """
