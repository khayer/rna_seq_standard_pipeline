from snakemake.utils import validate
import pandas as pd


validate(config, "../schemas/config.schema.yaml")

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
print(samples.head())
validate(samples, "../schemas/samples.schema.yaml")

samples = pd.read_csv(config["samples"], sep=",", dtype = str).set_index("Run", drop=False)

print(samples.head())
validate(samples, "../schemas/samples.schema.yaml")


def get_mem_gb(wildcards, attempt):
    return attempt * 60

def get_raw_reads(wildcards):
    run_ids = samples[samples["Run"] == wildcards.sample]["Run"].tolist()
    out = []
    if config["single_end"]:
        for r in run_ids:
            out.extend ( 
                expand ( [
                    "reads/{sample}.fastq.gz" 
                ], sample = r
                )
            )
    else:
        for r in run_ids:
            out.extend ( 
                expand ( [
                    "reads/{sample}_1.fastq.gz" , "reads/{sample}_2.fastq.gz"
                ], sample = r
                )
            )
    return out

def get_trimmed_reads(wildcards):
    run_ids = samples[samples["sample_name"] == wildcards.sample_name]["Run"].tolist()
    out = []
    if config["single_end"]:
        for r in run_ids:
            out.extend ( 
                expand ( [
                    "results/trimmed/{sample}_trim.fastq.gz" 
                ], sample = r
                )
            )
    else:
        for r in run_ids:
            out.extend ( 
                expand ( [
                    "results/trimmed/{sample}_trim_1.fastq.gz", "results/trimmed/{sample}_trim_2.fastq.gz"
                ], sample = r
                )
            )
    return out

def get_trimmed_reads2(wildcards):
    run_ids = samples[samples["Run"] == wildcards.sample]["Run"].tolist()
    out = []
    if config["single_end"]:
        for r in run_ids:
            out.extend ( 
                expand ( [
                    "results/trimmed/{sample}_trim.fastq.gz" 
                ], sample = r
                )
            )
    else:
        for r in run_ids:
            out.extend ( 
                expand ( [
                    "results/trimmed/{sample}_trim_1.fastq.gz", "results/trimmed/{sample}_trim_2.fastq.gz"
                ], sample = r
                )
            )
    return out


def get_sj_files(wildcards):
    run_ids = samples["sample_name"].tolist()
    out = []
    for r in run_ids:
        out.extend ( 
            expand ( [
                "results/splicing/majiq/majiq_{sample}/build_{sample}/{sample}.sj"
            ], sample = r
            )
        )
    return out

def get_all_star_junctions_files(wildcards):
    #print(wildcards)
    run_ids = samples["sample_name"].tolist()
    out = []
    for r in run_ids:
        out.extend ( 
            expand ( [
                "results/mapped/{sample}_SJ.out.tab"
            ], sample = r
            )
        )
    return out

def get_all_gene_tpm_files(wildcards):
    #print(wildcards)
    run_ids = samples["sample_name"].tolist()
    out = []
    for r in run_ids:
        out.extend ( 
            expand ( [
                "results/mapped/{sample}_Aligned.sortedByCoord.out_genes.out",
            ], sample = r
            )
        )
    return out

def get_all_gene_tpm_files_numbered_chr(wildcards):
    #print(wildcards)
    run_ids = samples["sample_name"].tolist()
    out = []
    for r in run_ids:
        out.extend ( 
            expand ( [
                "results/mapped/{sample}_numbered_chr_genes.out",
            ], sample = r
            )
        )
    return out


def get_all_read_counts(wildcards):
    #print(wildcards)
    run_ids = samples["sample_name"].tolist()
    out = []
    for r in run_ids:
        out.extend ( 
            expand ( [
                "results/mapped/{sample}_numbered_chr_number_of_reads.txt"
            ], sample = r
            )
        )
    return out

def get_all_gene_tpm_files_downsampled(wildcards):
    #print(wildcards)
    run_ids = samples["sample_name"].tolist()
    out = []
    for r in run_ids:
        out.extend ( 
            expand ( [
                "results/mapped_down/{sample}_numbered_chr_down_genes.out",
            ], sample = r
            )
        )
    return out



def all_input(wildcards):

    wanted_input = []

    ## QC with fastQC and multiQC
    #wanted_input.extend([
    #    "results/qc/multiqc/multiqc.html"
    #])

    # trimming
    for sample in samples.index:
        if config["single_end"]:
            wanted_input.extend(
                expand (
                        [
                            "results/trimmed/{sample}_trim.fastq.gz",
                            "results/fastqc/{sample}_trim_fastqc.zip"
                        ],
                        sample = sample
                    )
                )
        else:
            for read in ["1","2"]:
                wanted_input.extend(
                        expand (
                            [
                                "results/trimmed/{sample}_trim_{read}.fastq.gz",
                                "results/fastqc/{sample}_trim_{read}_fastqc.zip"
                            ],
                            sample = sample, read = read
                        )
                    )

    for sample in samples.sample_name:
        if not config["stranded"]:
            wanted_input.extend(
                expand (
                        [
                            "results/coverage/{sample}_CPM.bw"
                        ],
                        sample = sample
                    )
                )
        else:
            wanted_input.extend(
                    expand (
                        [
                            "results/coverage/{sample}_fwd_CPM.bw", 
                            "results/coverage/{sample}_rev_CPM.bw"
                        ],
                        sample = sample
                    )
                )

    for sn in samples.sample_name:
        wanted_input.extend(
                    expand (
                        [
                            
                            "results/quant/salmon_quant_{sample_name}/quant.sf",
                            #"results/splicing/majiq/majiq_{sample_name}/build_{sample_name}/{sample_name}.sj",
                            #"results/mapped/{sample_name}_Aligned.sortedByCoord.out_junc.bed",
                            "results/mapped/{sample_name}_Aligned.sortedByCoord.out_genes.ent",
                            "results/mapped/{sample_name}_selected_genes.bam",
                            "results/mapped/{sample_name}_numbered_chr.bam",
                            "results/mapped_down/{sample_name}_numbered_chr_down.bam",
                            "results/mapped/{sample_name}_idxstats.txt", 
                            "results/mapped/{sample_name}_flagstat.txt", 
                            "results/mapped/{sample_name}_samtools_stats.txt",
                            "results/mapped_down/{sample_name}_numbered_chr_down_intron_htscounts.txt"
                            
                        ],
                        sample_name = sn
                    )
                )
    ## get merged files
    wanted_input.extend(["results/quant/all_star_junctions.csv"])
    wanted_input.extend(["results/quant/all_gene_tmps.csv"])
    wanted_input.extend(["results/quant/all_gene_tmps_numbered_chr.csv"])
    wanted_input.extend(["results/quant/all_gene_tmps_numbered_chr_downsampled.csv"])
    return wanted_input
    



### run multiqc at the very end
rule run_multiqc:
    input: all_input
    output: "results/multiqc_report.html"
    log:    "00log/run_multiqc.log"
    conda: "../envs/python_tools.yaml"
    resources: 
        cpu = 2,
        mem = "10",
        time = "34:00:00"
    message: "run_multiqc {input[0]}: {resources.cpu} threads / {resources.mem}"
    shell:
        """
        cd results
        multiqc -f . 
        """