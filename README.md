---- UNDER CONSTRUCTION ----


# rna_seq_standard_pipeline

Snakemake workflow to analyze bulk RNA-seq. This is largerly inspired by https://github.com/snakemake-workflows.

## Samples config 

The sample spreadsheet specified in config.yaml has to be formatted like:

![example spreadsheet](misc/example_sample_table.png)

The config file should hold the following fields (this example is for hg38, but can be replaced for any other version/ species):

```bash
(snakemake) [hayerk@reslnvvhpc041 endpoints]$ head config/config.yaml
samples: config/SraRunTable.csv
star_index: /path/to/GENCODE_GRCh38.p13
salmon_index: /path/to/gencode.v30.transcripts_salmon
blacklist: /path/to/hg38-blacklist.v2.bed
gtf: /path/to/gencode.v32.primary_assembly.annotation.gtf
gff3: /path/to/gencode.v32.primary_assembly.annotation.gff3
tmp_dir: /scr1/users/hayerk/tmp
single_end: False
stranded: False
```

STAR is set to version 

Prepare STAR index example (make sure to provide the gtf file!):

    STAR --runMode genomeGenerate --genomeDir hg38_ad5_myco_dros_2.7.9a/ --genomeFastaFiles Adenovirus-Ad5.fasta GRCh38.p13/GRCh38.primary_assembly.genome.fa Mycoplasma/myco.fa dm6_fixed.fa --runThreadN 8 --sjdbGTFfile /mnt/isilon/thomas-tikhonenko_lab/data/index/GRCh38.p13_kat/gencode.v32.primary_assembly.annotation.gtf

## Folder structure

```bash
├── reads/ -> fastqfiles
├── config/
	├──samples.csv
	├── config.yaml
```

 ### Slurm (or other cluter) profile

    (snakemake) [hayerk@reslnvvhpc041 endpoints]$ head ~/.config/snakemake/slurm/config.yaml
    jobs: 500
    cluster: "sbatch -t {resources.time} --mem={resources.mem}G -c {resources.cpu} -J sm_{rule} -o logs_slurm/{rule}_{wildcards}.o -e logs_slurm/{rule}_{wildcards}.e" #--mail-type=FAIL --mail-user=hayerk@chop.edu"
    default-resources: [ mem=10, time=60, cpu=1]
    #resources: [cpus=30, mem_mb=500000]
	
## Testing 
    
    snakemake --profile slurm -s ~/path/to/rna_seq_standard_pipeline/workflow/Snakefile -p --use-conda --configfile config/config.yaml -n

## Submit to cluster

    snakemake --profile slurm -s ~/path/to/rna_seq_standard_pipeline/workflow/Snakefile -p --use-conda --configfile config/config.yaml


## Make this rulegraph

    snakemake -s workflow/Snakefile --configfile test_dir/config/config.yaml -d test_dir --rulegraph | dot -Tsvg > misc/dag.svg

![example dag](misc/dag.svg)


## How Kat runs it

    snakemake --cluster-cancel scancel --profile slurm -s ~/data/tools/rna_seq_standard_pipeline/workflow/Snakefile -p --use-conda --configfile config/config.yaml --rerun-triggers mtime --latency-wait 50 -k -n

## TODO

- Fix strand information (especially in regtools!)

