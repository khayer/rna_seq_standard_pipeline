rule trimming:
    input: "reads/{sample}_1.fastq.gz" , "reads/{sample}_2.fastq.gz"
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
    input: get_trimmed_reads
    output: "results/mapped/{sample_name}_Aligned.sortedByCoord.out.bam"
    log:    "00log/Star_align_{sample_name}.log"
    conda: "../envs/bioinf_tools.yaml"
    resources: 
        cpu = 10,
        mem = "60G",
        time = "44:00:00"
    params: 
        options = "--outFileNamePrefix {sample_name}_ --twopassMode Basic --outSAMtype BAM SortedByCoordinate --alignSJoverhangMin 8 --outSAMattributes All --outReadsUnmapped Fastx --readFilesCommand zcat --quantMode GeneCounts",
        star_index = config["star_index"],
        tmp_dir = config["tmp_dir"],
        sample_name = "{sample_name}"
    message: "aligning {input}: {resources.cpu} threads / {resources.mem}"
    shell:
        """
       	STAR {params.options} --genomeDir {params.star_index} --runThreadN {resources.cpu} --readFilesIn {input} --outTmpDir {params.tmp_dir}_{params.sample_name}
        samtools index {output}
        """

#rule salmon:
#    input: get_trimmed_reads
#    output: "results/mapped/{sample_name}_Aligned.sortedByCoord.out.bam"
#    log:    "00log/Star_align_{sample_name}.log"
#    conda: "../envs/bwa.yaml"
#    resources: 
#        cpu = 10,
#        mem = "40G",
#        time = "44:00:00"
#    params: 
#        options = "--twopassMode Basic --outSAMtype BAM SortedByCoordinate --alignSJoverhangMin 8 --outSAMattributes All --outReadsUnmapped Fastx --outFileNamePrefix {sample_name}_ --readFilesCommand zcat --quantMode GeneCounts",
#        star_index = config["star_index"],
#        tmp_dir = config["tmp_dir"],
#        sample_name = "{sample_name}"
#    message: "aligning {input}: {resources.cpu} threads / {resources.mem}"
#    shell:
#        """
#        STAR --genomeDir {params.star_index} --runThreadN {resources.cpu} --readFilesIn {input} --outTmpDir {params.tmp_dir}_{params.sample_name}
#        """

