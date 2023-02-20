import os
import sys
import pandas as pd
#from gtfparse import read_gtf


my_log = snakemake.log[0]
#my_log = "log_test.txt"

log_file = open(my_log, 'w')

def process_file(f, path):
    df =  pd.read_csv(path + f, sep = "\t", header =None, index_col = False, names = ["read_count"])
    #n = f.split("/")[1]
    n = f.replace("_numbered_chr_number_of_reads.txt","")
    #
    df["sample_name"] = n
    log_file.write(n + "\n")
    #print(df.head())
    #df.rename(columns={'TPM': 'TPM_' + n, 'Reads': 'Reads_' + n}, inplace=True)
    #df = df.iloc[:,4:6]
    #log_file.write(df.columns)
    
    return df


dir_mapped = snakemake.params.dir_mapped
#dir_mapped = sys.argv[1]
#print(dir_mapped)
log_file.write(dir_mapped + "\n")

outfile = snakemake.output[0]
#outfile = "results/downsample_spreadsheet.csv"
log_file.write(outfile + "\n")

sample_sheet = snakemake.params.sample_sheet
#sample_sheet = sys.argv[2]
sample_anno = pd.read_csv(sample_sheet)
print(sample_anno)


all_read_count_files = snakemake.input
#all_read_count_files = sys.argv[3:]
all_read_count_files = [sub.replace('results/mapped/', '') for sub in all_read_count_files]

print(all_read_count_files)
frames = [ process_file(f,dir_mapped) for f in all_read_count_files ]
result = pd.concat(frames, axis= 0, join='outer')

print(result.head())

final_sheet = sample_anno.merge(result)


final_sheet['min_by_background'] = final_sheet.groupby('background')['read_count'].transform('min')
final_sheet['ratio'] = final_sheet.groupby('background')['read_count'].transform(lambda x: x.min()/x)
print(final_sheet)

final_sheet.to_csv(outfile + "everything.csv" ,float_format = '{:.3f}'.format)
just_name_and_ratio = final_sheet[["sample_name","ratio"]]
just_name_and_ratio.to_csv(outfile  , sep = "\t", index = False,float_format = '{:.5f}'.format)
log_file.close()
