from snakemake.utils import validate
import pandas as pd

samples = pd.read_csv(config["samples"], sep=",", dtype = str).set_index("Run", drop=False)

def get_raw_reads(wildcards):
    run_ids = samples[samples["Run"] == wildcards.sample]["Run"].tolist()
    out = []
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
    for r in run_ids:
        out.extend ( 
            expand ( [
                "results/trimmed/{sample}_trim_1.fastq.gz", "results/trimmed/{sample}_trim_2.fastq.gz"
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
        for read in ["1","2"]:
            wanted_input.extend(
                    expand (
                        [
                            "results/trimmed/{sample}_trim_{read}.fastq.gz"
                        ],
                        sample = sample, read = read
                    )
                )

    for sn in samples.sample_name:
        wanted_input.extend(
                    expand (
                        [
                            "results/mapped/{sample_name}_Aligned.sortedByCoord.out.bam",
                            "results/quant/salmon_quant_{sample_name}/quant.sf",
                            "results/splicing/majiq/majiq_{sample_name}/build_{sample_name}/{sample_name}.sj"#,
                            #"results/coverage/{sample_name}_fwd_CPM.bw", 
                            #"results/coverage/{sample_name}_rev_CPM.bw"
                        ],
                        sample_name = sn
                    )
                )

    return wanted_input
    