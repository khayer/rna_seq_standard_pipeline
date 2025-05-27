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
        regtools junctions extract -s RF -o {output[0]} {input[0]}
        """

rule run_htseq:
    input: "results/mapped/{sample_name}_Aligned.sortedByCoord.out.bam", "results/mapped/{sample_name}_Aligned.sortedByCoord.out.bam.bai"
    output: "results/mapped/{sample_name}_htscounts.txt"
    log:    "00log/run_htseq_{sample_name}.log"
    conda: "../envs/htseq_env.yaml"
    resources: 
        cpu = 14,
        mem = "20",
        time = "24:00:00"
    params: 
        gtf_anno = config["gtf"]
    message: "run_htseq {input}: {resources.cpu} threads / {resources.mem}"
    shell:
        """
        htseq-count -n {resources.cpu} -f bam -t exon -r pos -s reverse -m intersection-strict {input[0]} {params.gtf_anno} > {output[0]}
        """

rule run_htseq_downsampled:
    input: "results/mapped_down/{sample_name}_numbered_chr_down.bam", "results/mapped_down/{sample_name}_numbered_chr_down.bam.bai"
    output: "results/mapped_down/{sample_name}_htscounts_down.txt"
    log:    "00log/run_htseq_{sample_name}.log"
    conda: "../envs/htseq_env.yaml"
    resources: 
        cpu = 14,
        mem = "20",
        time = "24:00:00"
    params: 
        gtf_anno = config["gtf"]
    message: "run_htseq_downsampled {input}: {resources.cpu} threads / {resources.mem}"
    shell:
        """
        htseq-count -n {resources.cpu} -f bam -t exon -r pos -s reverse -m intersection-strict {input[0]} {params.gtf_anno} > {output[0]}
        """

### calculate TPMs
if config["single_end"]:
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
            TPMCalculator -g {params.gtf_anno} -b {params.in_file} -q 200 -e 
            """

    rule run_TPMCalculator_numbered_chr:
        input: "results/mapped/{sample_name}_numbered_chr.bam", "results/mapped/{sample_name}_numbered_chr.bam.bai"
        output: "results/mapped/{sample_name}_numbered_chr_genes.ent","results/mapped/{sample_name}_numbered_chr_genes.out"
        log:    "00log/run_TPMCalculator-numbered_chr_{sample_name}.log"
        conda: "../envs/bioinf_tools.yaml"
        resources: 
            cpu = 2,
            mem = "10",
            time = "34:00:00"
        #wildcard_constraints:
        #    sample_name="\w+"
        params: 
            gtf_anno = config["gtf"],
            in_file = "{sample_name}_numbered_chr.bam"
        message: "run_TPMCalculator_numbered_chr {input}: {resources.cpu} threads / {resources.mem}"
        shell:
            """
            cd results/mapped/
            TPMCalculator -g {params.gtf_anno} -b {params.in_file} -q 200 -e 
            """

    rule run_TPMCalculator_downsampled:
        input: "results/mapped_down/{sample_name}_numbered_chr_down.bam", "results/mapped_down/{sample_name}_numbered_chr_down.bam.bai"
        output: "results/mapped_down/{sample_name}_numbered_chr_down_genes.ent","results/mapped_down/{sample_name}_numbered_chr_down_genes.out"
        log:    "00log/run_TPMCalculator_downsampled_{sample_name}.log"
        conda: "../envs/bioinf_tools.yaml"
        resources: 
            cpu = 2,
            mem = "10",
            time = "34:00:00"
        params: 
            gtf_anno = config["gtf"],
            in_file = "{sample_name}_numbered_chr_down.bam"
        message: "rrun_TPMCalculator_downsampled {input}: {resources.cpu} threads / {resources.mem}"
        shell:
            """
            cd results/mapped_down/
            TPMCalculator -g {params.gtf_anno} -b {params.in_file} -q 200 -e 
            """

else:
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

    rule run_TPMCalculator_numbered_chr:
        input: "results/mapped/{sample_name}_numbered_chr.bam", "results/mapped/{sample_name}_numbered_chr.bam.bai"
        output: "results/mapped/{sample_name}_numbered_chr_genes.ent","results/mapped/{sample_name}_numbered_chr_genes.out"
        log:    "00log/run_TPMCalculator_numbered_chr_{sample_name}.log"
        conda: "../envs/bioinf_tools.yaml"
        resources: 
            cpu = 2,
            mem = "10",
            time = "34:00:00"
        params: 
            gtf_anno = config["gtf"],
            in_file = "{sample_name}_numbered_chr.bam"
        message: "run_TPMCalculator_numbered_chr {input}: {resources.cpu} threads / {resources.mem}"
        shell:
            """
            cd results/mapped/
            TPMCalculator -g {params.gtf_anno} -b {params.in_file} -p -q 200 -e 
            """

    rule run_TPMCalculator_downsampled:
        input: "results/mapped_down/{sample_name}_numbered_chr_down.bam", "results/mapped_down/{sample_name}_numbered_chr_down.bam.bai"
        output: "results/mapped_down/{sample_name}_numbered_chr_down_genes.ent","results/mapped_down/{sample_name}_numbered_chr_down_genes.out"
        log:    "00log/run_TPMCalculator_downsampled_{sample_name}.log"
        conda: "../envs/bioinf_tools.yaml"
        resources:
            cpu = 2,
            mem = "10",
            time = "34:00:00"
        params:
            gtf_anno = config["gtf"],
            in_file = "{sample_name}_numbered_chr_down.bam"
        message: "run_TPMCalculator_downsampled {input}: {resources.cpu} threads / {resources.mem}"
        shell:
            """
            cd results/mapped_down/
            TPMCalculator -g {params.gtf_anno} -b {params.in_file} -p -q 200 -e 
            """


rule run_htseq_downsampled_introns:
    input: "results/mapped_down/{sample_name}_numbered_chr_down.bam", "results/mapped_down/{sample_name}_numbered_chr_down.bam.bai"
    output: "results/mapped_down/{sample_name}_numbered_chr_down_intron_htscounts.txt"
    log:    "00log/run_TPMCalculator_downsampled_{sample_name}.log"
    conda: "../envs/htseq_env.yaml"
    resources: 
        cpu = 8,
        mem = "20",
        time = "24:00:00"
    params: 
        gtf_anno_introns = config["gtf_introns"]
    message: "run_htseq_downsampled_introns {input}: {resources.cpu} threads / {resources.mem}"
    shell:
        """
        htseq-count -n {resources.cpu} -f bam -t exon -r pos -s no -m intersection-strict {input[0]} {params.gtf_anno_introns} > {output[0]}
        """


rule run_merge_junctions_STAR:
    input: get_all_star_junctions_files
    output: "results/quant/all_star_junctions.csv"
    log:    "00log/run_merge_junctions_STAR.log"
    conda: "../envs/bioinf_tools.yaml"
    resources: 
        cpu = 2,
        mem = "40",
        time = "44:00:00"
    params: "results/mapped/"
    message: "run_merge_junctions_STAR {params}: {resources.cpu} threads / {resources.mem}"
    shell:
        "python3 {workflow.basedir}/scripts/merge_STAR_junctions.py {params} {output} {input}"

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

rule run_merge_gene_tpms_numbered_chr:
    input: get_all_gene_tpm_files_numbered_chr
    output: "results/quant/all_gene_tmps_numbered_chr.csv"
    log:    "00log/run_merge_gene_tpms_numbered_chr.log"
    conda: "../envs/python_tools.yaml"
    resources: 
        cpu = 2,
        mem = "30",
        time = "24:00:00"
    params: "results/mapped/", config["gtf"]
    message: "run_merge_gene_tpms_numbered_chr {params}: {resources.cpu} threads / {resources.mem}"
    script:
        "../scripts/merge_gene_TPMs.py"

rule run_merge_gene_tpms_numbered_downsampled:
    input: get_all_gene_tpm_files_downsampled
    output: "results/quant/all_gene_tmps_numbered_chr_downsampled.csv"
    log:    "00log/run_merge_gene_tpms_numbered_downsampled.log"
    conda: "../envs/python_tools.yaml"
    resources: 
        cpu = 2,
        mem = "30",
        time = "24:00:00"
    params: "results/mapped_down/", config["gtf"]
    message: "run_merge_gene_tpms_numbered_downsampled {params}: {resources.cpu} threads / {resources.mem}"
    script:
        "../scripts/merge_gene_TPMs.py"

rule run_merge_gene_htseq_counts:
    input: get_all_gene_htseq_files
    output: "results/quant/all_gene_htseq_counts.csv"
    log:    "00log/run_merge_gene_htseq_counts.log"
    conda: "../envs/python_tools.yaml"
    resources: 
        cpu = 2,
        mem = "30",
        time = "24:00:00"
    params: "results/mapped/", config["gtf"]
    message: "run_merge_gene_htseq_counts {params}: {resources.cpu} threads / {resources.mem}"
    script:
        "../scripts/merge_htseq.py"


rule run_merge_gene_htseq_counts_down:
    input: get_all_gene_htseq_files_down
    output: "results/quant/all_gene_htseq_counts_numbered_chr_downsampled.csv"
    log:    "00log/run_merge_gene_htseq_counts.log"
    conda: "../envs/python_tools.yaml"
    resources: 
        cpu = 2,
        mem = "30",
        time = "24:00:00"
    params: "results/mapped_down/", config["gtf"]
    message: "run_merge_gene_htseq_counts_down {params}: {resources.cpu} threads / {resources.mem}"
    script:
        "../scripts/merge_htseq.py"
