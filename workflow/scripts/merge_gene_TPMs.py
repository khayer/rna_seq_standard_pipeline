import os
import sys
import pandas as pd
from gtfparse import read_gtf


log_file = open(snakemake.log[0], 'w')

def process_file(f, path):
    df =  pd.read_csv(path + f, sep = "\t",index_col = 0)
    #n = f.split("/")[1]
    n = f.replace("_Aligned.sortedByCoord.out_genes.out","")
    #
    log_file.write(n + "\n") 
    df.rename(columns={'TPM': 'TPM_' + n, 'Reads': 'Reads_' + n}, inplace=True)
    df = df.iloc[:,4:6]
    #log_file.write(df.columns)
    
    return df


dir_TPM = snakemake.params[0]
log_file.write(dir_TPM + "\n")

outfile = snakemake.output[0]
log_file.write(outfile + "\n")
all_TPM_files = [sub.replace('results/mapped/', '') for sub in snakemake.input]


anno = read_gtf(snakemake.params[1])

# filter DataFrame to gene entries on chrY
anno_genes = anno[anno["feature"] == "gene"]
anno_genes_chrY = anno_genes[anno_genes["seqname"] == "Y"]
#anno.head()
anno_genes = anno_genes[anno_genes.columns[0:11]]
anno_genes.index = anno_genes.gene_id


log_file.write("got passed reading annotation\n")

frames = [ process_file(f,dir_TPM) for f in all_TPM_files ]
result = pd.concat(frames, axis=1, join='outer')

log_file.write("got passed merging results\n")

# write to file 
#results_final = result
filter_col = [col for col in result if col.startswith('TPM')]
result[filter_col].join(anno_genes).to_csv(outfile,float_format = '{:.0f}'.format)

filter_col = [col for col in result if col.startswith('Reads_')]
result[filter_col].join(anno_genes).to_csv(outfile + "_Reads.csv",float_format = '{:.0f}'.format)

log_file.close()
