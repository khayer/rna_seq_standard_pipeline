from snakemake.utils import validate
import pandas as pd

samples = pd.read_csv(config["samples"], sep=",", dtype = str).set_index("Run", drop=False)


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

    return wanted_input
    