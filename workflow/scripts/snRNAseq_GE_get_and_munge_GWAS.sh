#!/bin/bash

# SCZ data wave 3 v2 file is here: https://doi.org/10.6084/m9.figshare.14672178
# Height data Wood et al + UK BioBank is here: https://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files

# Set variables
SCZ_LOC=https://figshare.com/ndownloader/files/28169757
HEIGHT_LOC=https://portals.broadinstitute.org/collaboration/giant/images/6/63/Meta-analysis_Wood_et_al%2BUKBiobank_2018.txt.gz
OUT_DIR=../results/GWAS_for_MAGMA

mkdir ${OUT_DIR}
cd ${OUT_DIR}

# Download GWAS
wget ${SCZ_LOC}
wget ${HEIGHT_LOC}

# Unpack GWAS
mv 28169757 PGC3_SCZ_wave3_public.v2.tsv.gz
gunzip Meta-analysis_Wood_et_al+UKBiobank_2018.txt.gz
gunzip PGC3_SCZ_wave3_public.v2.tsv.gz 

# Get GWAS ready for MAGMA
awk '{print $3"\t"$1"\t"$2"\t"$9"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$10}' Meta-analysis_Wood_et_al+UKBiobank_2018.txt |\
sed 's/POS/BP/g' > HEIGHT_hg19_magma_ready.tsv

awk '{print $2"\t"$1"\t"$3"\t"$11"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' PGC3_SCZ_wave3_public.v2.tsv |\
sed 's/POS/BP/g' > SCZ_hg19_magma_ready.tsv

rm PGC3_SCZ_wave3_public.v2.tsv Meta-analysis_Wood_et_al+UKBiobank_2018.txt
