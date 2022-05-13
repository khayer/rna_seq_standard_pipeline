### extract junctions
rule run_regtools:
    input: "results/mapped/{sample_name}_Aligned.sortedByCoord.out.bam", "results/mapped/{sample_name}_Aligned.sortedByCoord.out.bam.bai"
    output: "results/mapped/{sample_name}_Aligned.sortedByCoord.out_junc.bed"
    log: "00log/run_regtools_{sample_name}.log"
    conda: "../envs/bioinf_tools.yaml"
    resources: 
        cpu = 2,
        mem = "10",
        time = "24:00:00"
    message: "run_regtools {input}: {resources.cpu} threads / {resources.mem}"
    shell:
        """
        regtools junctions extract -s 0 -o {output[0]} {input[0]}
        """

### calculate TPMs
rule run_TPMCalculator:
    input: "results/mapped/{sample_name}_Aligned.sortedByCoord.out.bam", "results/mapped/{sample_name}_Aligned.sortedByCoord.out.bam.bai"
    output: "results/mapped/{sample_name}_Aligned.sortedByCoord.out_genes.ent","results/mapped/{sample_name}_Aligned.sortedByCoord.out_genes.out"
    log:    "00log/run_TPMCalculator_{sample_name}.log"
    conda: "../envs/bioinf_tools.yaml"
    resources: 
        cpu = 2,
        mem = "10",
        time = "34:00:00"
    params: 
        gtf_anno = config["gtf"],
        in_file = "{sample_name}_Aligned.sortedByCoord.out.bam"
    message: "run_TPMCalculator {input}: {resources.cpu} threads / {resources.mem}"
    shell:
        """
        cd results/mapped/
        TPMCalculator -g {params.gtf_anno} -b {params.in_file} -p -q 200 -e 
        """

rule run_merge_junctions_STAR:
    input: get_all_star_junctions_files
    output: "results/quant/all_star_junctions.csv"
    log:    "00log/run_merge_junctions_STAR.log"
    #conda: "../envs/bioinf_tools.yaml"
    resources: 
        cpu = 2,
        mem = "20",
        time = "24:00:00"
    params: "results/mapped/"
    message: "run_merge_junctions_STAR {params}: {resources.cpu} threads / {resources.mem}"
    script:
        "../scripts/merge_STAR_junctions.py"

rule run_merge_gene_tpms:
    input: get_all_gene_tpm_files
    output: "results/quant/all_gene_tmps.csv"
    log:    "00log/run_merge_gene_tpms.log"
    conda: "../envs/python_tools.yaml"
    resources: 
        cpu = 2,
        mem = "30",
        time = "24:00:00"
    params: "results/mapped/", config["gtf"]
    message: "run_merge_gene_tpms {params}: {resources.cpu} threads / {resources.mem}"
    script:
        "../scripts/merge_gene_TPMs.py"

