import os
import sys
import pandas as pd

log_file = open(snakemake.log[0], 'w')

#args = sys.argv
#args[1] should hold the directory with junction files
dir_junc = snakemake.params[0]
log_file.write(dir_junc + "\n")

#args[2] should name and directory of the output file
outfile = snakemake.output[0]
log_file.write(outfile + "\n")
#all_junc_files = os.listdir(dir_junc)
#all_junc_files = [name for name in all_junc_files if '_SJ.out.tab' in name]
#len(all_junc_files)
all_junc_files = snakemake.input
all_junc_files = [sub.replace('results/mapped/', '') for sub in snakemake.input]
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

log_file.write("ALL GOOD\n")

#junc_df.to_csv(outfile + "_kat_version.csv")

## CHANGE FORMAT
#(snakemake) [hayerk@reslnvvhpc040 cluster_training]$ head results/quant/all_star_junctions.csv
#,chr,start,stop,unique,norm,sample_name
#5,chr1,14830,14969,110,3.3788617548542415,SLN2420
#10,chr1,15039,15795,30,0.9215077513238841,SLN2420
#11,chr1,15174,185616,6,0.1843015502647768,SLN2420
#14,chr1,15948,16606,6,0.1843015502647768,SLN2420
#17,chr1,16766,16857,7,0.21501847530890628,SLN2420
#24,chr1,17526,188049,10,0.3071692504412947,SLN2420

pivot_df = junc_df.pivot_table(index = ["chr","start", "stop"], columns="sample_name", values = "norm").fillna(0)
log_file.write("Pivot done\n")

pivot_df.to_csv(outfile)
log_file.close()