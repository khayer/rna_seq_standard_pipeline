import os
import pandas as pd
import sys


#args = sys.argv
#args[1] should hold the directory with junction files
dir_junc = snakemake.params[0]
print(dir_junc)

#args[2] should name and directory of the output file
outfile = snakemake.output[0]
print(outfile)
all_junc_files = os.listdir(dir_junc)
all_junc_files = [name for name in all_junc_files if '_SJ.out.tab' in name]
#len(all_junc_files)


junctions = pd.read_csv(dir_junc + all_junc_files[0], sep = "\t", header = None,
	usecols = [0,1,2,3,5,6,7,8],
	names = ["chr","start","stop","strand","annotated","unique","multi","max"])
#junctions.head()
all_reads = junctions['unique'].sum()
junctions["norm"] = junctions["unique"]/(all_reads / 1000000)

junc_df = junctions.loc[:,("chr","start","stop","unique","norm")]
junc_df["sample_name"] = "test"
junc_df.head()

junc_df = pd.DataFrame(columns=junc_df.columns)
#for f in all_voila_files[0:5]:

for f in all_junc_files:

    #print(f)
    junctions.head()
    c_junc_df = pd.read_csv(dir_junc+f,
                            sep='\t',
                            header=None,
                            usecols = [0,1,2,3,5,6,7,8],
                            names = ["chr","start","stop","strand","annotated","unique","multi","max"]
                           )
    all_reads = c_junc_df['unique'].sum()
    c_junc_df = c_junc_df[c_junc_df['unique'] > 5]
    c_junc_df['norm'] = c_junc_df['unique']/(all_reads/1000000)
    
    c_junc_df["sample_name"] = f.split("_SJ.out.tab")[0]
    
   
    c_junc_df = c_junc_df.loc[:,("chr","start","stop","unique","norm","sample_name")]
    junc_df  = pd.concat([junc_df, c_junc_df])
    

junc_df.SRA.unique().shape
    

junc_df.shape

junc_df.to_csv(outfile)


