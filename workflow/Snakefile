# The main entry point of workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.
# taken from https://github.com/snakemake-workflows
import sys
include: "rules/common.smk"
#include: "rules/ref.smk"
#include: "rules/qc.smk"
#include: "rules/cutadapt.smk"
include: "rules/mapping.smk"
include: "rules/splicing.smk"

rule all:
    input: all_input

if not os.path.exists("00log"):
    os.makedirs("00log")
if not os.path.exists("logs_slurm"):
    os.makedirs("logs_slurm")
if not os.path.exists("results"):
    os.makedirs("results")
if not os.path.exists("results/coverage"):
    os.makedirs("results/coverage")
