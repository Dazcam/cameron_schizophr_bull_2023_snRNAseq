#!/usr/bin/env

mkdir -p ../resources/raw_data/shi_et_al_2021

# Metadata
wget https://www.science.org/doi/suppl/10.1126/science.abj6641/suppl_file/science.abj6641_tables_s2_to_s9.zip -P ../resources/raw_data/shi_et_al_2021
gunzip ../resources/raw_data/shi_et_al_2021/science.abj6641_tables_s2_to_s9.zip

# GeX matrix
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135827/suppl/GSE135827%5FGE%5Fmat%5Fraw%5Fcount%5Fwith%5Fweek%5Finfo%2Etxt%2Egz

