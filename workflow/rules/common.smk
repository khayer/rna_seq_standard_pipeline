from snakemake.utils import validate
import pandas as pd

samples = pd.read_csv(config["samples"], sep=",", dtype = str).set_index("Run", drop=False)
samples.index.names = ["sample_id"]


def get_mapped_reads(wildcards):
    run_ids = samples[samples["kind"] == wildcards.group]["Run"].tolist()
    out = []
    for r in run_ids:
        out.extend ( 
            expand ( [
                "results/mapped/{sample}.unsorted_{read}.bam"
            ], sample = r, read = wildcards.read
            )
        )
    return out

def get_mapped_reads_formatted(wildcards):
    run_ids = samples[samples["kind"] == wildcards.group]["Run"].tolist()
    out = []
    for r in run_ids:
        out.extend ( 
            expand ( [
                "results/mapped/{sample}.unsorted_{read}.bam"
            ], sample = r, read = wildcards.read
            )
        )
    out = " ".join(["-INPUT " + sub for sub in out])
    return out



def all_input(wildcards):

    wanted_input = []

    ## QC with fastQC and multiQC
    #wanted_input.extend([
    #    "results/qc/multiqc/multiqc.html"
    #])

    # mapping
    for sample in samples.index:
        for read in ["1","2"]:
            wanted_input.extend(
                    expand (
                        [
                            "results/mapped/{sample}.unsorted_{read}.bam"
                        ],
                        sample = sample, read = read
                    )
                )

    for kind in samples["kind"].unique().tolist():
        for read in ["1","2"]:
            wanted_input.extend(
                    expand (
                        [
                            "results/merged/{group}.unsorted_{read}.bam"
                        ],
                        group = kind, read = read
                    )
                )
    for kind in samples["kind"].unique().tolist():
        wanted_input.extend(
            expand (
                [ "results/hic_matrix/raw_{group}_bs_1000_master.h5"
                ], group = kind 
                )
    )
    return wanted_input
    