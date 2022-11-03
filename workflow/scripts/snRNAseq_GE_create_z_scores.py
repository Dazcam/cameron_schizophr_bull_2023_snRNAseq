import pandas as pd

genes_scz = pd.read_csv(snakemake.input[0], delim_whitespace = True)
genes_height = pd.read_csv(snakemake.input[1], delim_whitespace = True)
mag = pd.read_csv(snakemake.input[2], delimiter="\t")
mag.columns = ['GENE', 'chr', 'start', 'stop', 'strand', 'SYMBOL']
mag_filt = mag[['GENE','SYMBOL']]
scz_merge = mag_filt.merge(genes_scz, how="left", on=['GENE', 'GENE'])
scz_merge = scz_merge[['SYMBOL', 'ZSTAT', 'GENE']]
scz_merge.columns = ['SYMBOL','SCZ', 'GENE']
height_merge = scz_merge.merge(genes_height, how="left", on=['GENE', 'GENE'])
height_merge = height_merge[['SYMBOL', 'SCZ', 'ZSTAT']]
height_merge.columns = ['GENE', 'SCZ', 'HEIGHT']

height_merge.to_csv(snakemake.output[0], header = True, index = False, sep = '\t')
