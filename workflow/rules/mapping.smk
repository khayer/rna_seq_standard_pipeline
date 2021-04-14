rule align:
    input: "reads/{sample}_{read}.fastq.gz"
    output: "results/mapped/{sample}.unsorted_{read}.bam"
    log:    "00log/bwa_align_{sample}_{read}.log"
    conda: "../envs/bwa.yaml"
    resources: 
        cpu = 10,
        mem = "40G",
        time = "44:00:00"
    params: 
        bwa = "-A1 -B4  -E50 -L0" 
    message: "aligning {input}: {resources.cpu} threads / {resources.mem}"
    shell:
        """
       	bwa mem {params.bwa} -t {resources.cpu} {config[idx_bwa]} {input} 2> {log} |  samtools view -Shb - >  {output}
        """




rule merge_bam_by_kind:
    input: get_mapped_reads
    output: "results/merged/{group}.unsorted_{read}.bam"
    log:    "00log/merge_bam_by_kind_{group}_{read}.log"
    conda: "../envs/picard.yaml"
    resources: 
        cpu = 10,
        mem = "20G",
        time = "44:00:00"
    params: 
        picard = "-SORT_ORDER queryname",
        tmp_dir = config["tmp_dir"],
        input = get_mapped_reads_formatted
    message: "merge_bam_by_kind {input}: {resources.cpu} threads / {resources.mem}"
    shell:
        """
        picard MergeSamFiles {params.picard} {params.input} -OUTPUT {output} -TMP_DIR {params.tmp_dir} 2> {log}  
        """
