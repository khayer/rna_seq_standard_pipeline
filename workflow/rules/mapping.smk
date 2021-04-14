rule trimming:
    input: "reads/{sample}_1.fastq.gz", "reads/{sample}_2.fastq.gz"
    output: "results/trimmed/{sample}_trim_1.fastq.gz", "results/trimmed/{sample}_trim_2.fastq.gz"
    log:    "00log/trim_{sample}.log"
    conda: "../envs/bioinf_tools.yaml"
    resources: 
        cpu = 10,
        mem = "20G",
        time = "44:00:00"
    params: 
        options = "ktrim=r k=23 mink=11 hdist=1 minlength=35 tpe tbo qtrim=r trimq=20 qin=33" 
    message: "trimming {input}: {resources.cpu} threads / {resources.mem}"
    shell:
        """
        bbduk.sh in={input[0]} in2={input[1]} ref=adapters {params.options} -Xmx10g threads={resources.cpu} out={output[0]} out2={output[1]} 2> {log} 
        """

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

