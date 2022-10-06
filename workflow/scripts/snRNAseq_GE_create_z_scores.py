import pandas as pd

genes = pd.read_csv(snakemake.input[0], delim_whitespace = True)
mag = pd.read_csv(snakemake.input[1], delimiter="\t")
mag.columns = ['GENE', 'chr', 'start', 'stop', 'strand', 'SYMBOL']
mag_filt = mag[['GENE','SYMBOL']]
left_merged = pd.merge(genes, mag_filt, how="left", on=['GENE', 'GENE'])
left_merged = left_merged[['SYMBOL','ZSTAT']]
left_merged.columns = ['GENE', 'SCZ']

left_merged.to_csv(snakemake.output[0], header = True, index = False, sep = '\t')
