rule trimming:
    input: "reads/{sample}_1.fastq.gz" , "reads/{sample}_2.fastq.gz"
    output: "results/trimmed/{sample}_trim_1.fastq.gz", "results/trimmed/{sample}_trim_2.fastq.gz"
    log:    "00log/trim_{sample}.log"
    conda: "../envs/bioinf_tools.yaml"
    resources: 
        cpu = 10,
        mem = "20",
        time = "44:00:00"
    params: 
        options = "ktrim=r k=23 mink=11 hdist=1 minlength=35 tpe tbo qtrim=r trimq=20 qin=33" 
    message: "trimming {input}: {resources.cpu} threads / {resources.mem}"
    shell:
        """
        bbduk.sh in={input[0]} in2={input[1]} ref=adapters {params.options} -Xmx10g threads={resources.cpu} out={output[0]} out2={output[1]} 2> {log} 
        """

rule trimming_single:
    input: "reads/{sample}.fastq.gz" 
    output: "results/trimmed/{sample}_trim.fastq.gz"
    log:    "00log/trim_{sample}.log"
    conda: "../envs/bioinf_tools.yaml"
    resources: 
        cpu = 10,
        mem = "20",
        time = "44:00:00"
    params: 
        options = "ktrim=r k=23 mink=11 hdist=1 minlength=35 tpe tbo qtrim=r trimq=20 qin=33" 
    message: "trimming {input}: {resources.cpu} threads / {resources.mem}"
    shell:
        """
        bbduk.sh in={input[0]} ref=adapters {params.options} -Xmx10g threads={resources.cpu} out={output[0]} 2> {log} 
        """


rule fastqc:
    input: get_trimmed_reads2
    output: "results/fastqc/{sample}_trim_1_fastqc.zip","results/fastqc/{sample}_trim_2_fastqc.zip"
    log:    "00log/fastqc_{sample}.log"
    conda: "../envs/bioinf_tools.yaml"
    resources: 
        cpu = 6,
        mem = "10",
        time = "24:00:00"
    params: 
        options = " "
    message: "fastqc {input}: {resources.cpu} threads / {resources.mem}"
    shell:
        """
        fastqc -t {resources.cpu} -o results/fastqc/ {input}
        """

rule fastqc_single:
    input: get_trimmed_reads2
    output: "results/fastqc/{sample}_trim_fastqc.zip"
    log:    "00log/fastqc_{sample}.log"
    conda: "../envs/bioinf_tools.yaml"
    resources: 
        cpu = 6,
        mem = "10",
        time = "24:00:00"
    params: 
        options = " "
    message: "fastqc {input}: {resources.cpu} threads / {resources.mem}"
    shell:
        """
        fastqc -t {resources.cpu} -o results/fastqc/ {input}
        """

def get_mem_gb_star(wildcards, attempt):
    return ((attempt * 60) - 10 ) * 1024 * 1024

rule align:
    input: get_trimmed_reads
    output: "results/mapped/{sample_name}_Aligned.sortedByCoord.out.bam", "results/mapped/{sample_name}_Aligned.sortedByCoord.out.bam.bai", "results/mapped/{sample_name}_SJ.out.tab"
    log:    "00log/Star_align_{sample_name}.log"
    conda: "../envs/bioinf_tools.yaml"
    resources: 
        cpu = 10,
        mem = get_mem_gb,
        time = "44:00:00",
        sort_mem = get_mem_gb_star
    params: 
        options = "--outFileNamePrefix results/mapped/{sample_name}_ --twopassMode Basic --outSAMtype BAM Unsorted --alignSJoverhangMin 8 --outSAMattributes All --outReadsUnmapped Fastx --readFilesCommand zcat --quantMode GeneCounts",
        star_index = config["star_index"],
        tmp_dir = config["tmp_dir"],
        sample_name = "{sample_name}"
        
    message: "aligning {input}: {resources.cpu} threads / {resources.mem}"
    shell:
        """
        if [ -d {params.tmp_dir}/STAR_{params.sample_name} ]; then rm -r {params.tmp_dir}/STAR_{params.sample_name}; fi
        #STAR {params.options} --limitBAMsortRAM {resources.sort_mem} --genomeDir {params.star_index} --runThreadN {resources.cpu} --readFilesIn {input} --outTmpDir {params.tmp_dir}/STAR_{params.sample_name}
        STAR {params.options} --genomeDir {params.star_index} --runThreadN {resources.cpu} --readFilesIn {input} --outTmpDir {params.tmp_dir}/STAR_{params.sample_name}
        samtools sort -T {params.tmp_dir}/samtools_{params.sample_name} -m 3G -@ {resources.cpu} results/mapped/{params.sample_name}_Aligned.out.bam > {output[0]}
        samtools index {output[0]} 
        rm results/mapped/{params.sample_name}_Aligned.out.bam
        rm -r {params.tmp_dir}/STAR_{params.sample_name}
        """


if config["single_end"]:
    rule salmon_single:
        input: get_trimmed_reads
        output: "results/quant/salmon_quant_{sample_name}/quant.sf"
        log:    "00log/Salmon_quant_{sample_name}.log"
        conda: "../envs/salmon.yaml"
        resources: 
            cpu = 10,
            mem = "40",
            time = "12:00:00"
        params: 
            salmon_index = config["salmon_index"],
            tmp_dir = config["tmp_dir"],
            sample_name = "{sample_name}"
        message: "salmon quant {input}: {resources.cpu} threads / {resources.mem}"
        shell:
            """
            salmon quant -i {params.salmon_index} -p {resources.cpu} -l A --validateMappings -r {input[0]} -o results/quant/salmon_quant_{params.sample_name}
            """

else:
    rule salmon:
        input: get_trimmed_reads
        output: "results/quant/salmon_quant_{sample_name}/quant.sf"
        log:    "00log/Salmon_quant_{sample_name}.log"
        conda: "../envs/salmon.yaml"
        resources: 
            cpu = 10,
            mem = "40",
            time = "12:00:00"
        params: 
            salmon_index = config["salmon_index"],
            tmp_dir = config["tmp_dir"],
            sample_name = "{sample_name}"
        message: "salmon quant {input}: {resources.cpu} threads / {resources.mem}"
        shell:
            """
            salmon quant -i {params.salmon_index} -p {resources.cpu} -l A --validateMappings -1 {input[0]} -2 {input[1]} -o results/quant/salmon_quant_{params.sample_name}
            """


if config["single_end"]:
    rule bamCoverage_CPM_single:
        input: "results/mapped/{sample_name}_Aligned.sortedByCoord.out.bam", "results/mapped/{sample_name}_Aligned.sortedByCoord.out.bam.bai"
        output: "results/coverage/{sample_name}_CPM.bw"
        log:  "00log/{sample_name}.bamCoverage_CPM_single"
        conda: "../envs/deeptools.yaml"
        resources:
            cpu = 4,
            mem = "10",
            time = "12:00:00"
        params:
            blacklist = "--blackListFileName " + config["blacklist"]
        message: "bamCoverage_CPM_single {input}: {resources.cpu} threads" #"/ {params.mem}"
        shell: 
            """
            bamCoverage -bs 1 -b {input[0]} -o {output[0]} -p 4 --normalizeUsing CPM {params.blacklist} --exactScaling --ignoreDuplicates --minMappingQuality 255
            """


else:
    rule bamCoverage_CPM:
        input: "results/mapped/{sample_name}_Aligned.sortedByCoord.out.bam", "results/mapped/{sample_name}_Aligned.sortedByCoord.out.bam.bai"
        output: "results/coverage/{sample_name}_fwd_CPM.bw", "results/coverage/{sample_name}_rev_CPM.bw"
        log:  "00log/{sample_name}.bamCoverage"
        conda: "../envs/deeptools.yaml"
        resources:
            cpu = 4,
            mem = "10",
            time = "12:00:00"
        params:
            blacklist = "--blackListFileName " + config["blacklist"]
        message: "CPM_bamCoverage_fwd_rev {input}: {resources.cpu} threads" #"/ {params.mem}"
        shell: 
            """
            bamCoverage -bs 1 -b {input[0]} -o {output[0]} --filterRNAstrand forward -p 4 --normalizeUsing CPM {params.blacklist} --exactScaling --ignoreDuplicates --minMappingQuality 255
            bamCoverage -bs 1 -b {input[0]} -o {output[1]} --filterRNAstrand reverse -p 4 --normalizeUsing CPM {params.blacklist} --exactScaling --ignoreDuplicates --minMappingQuality 255
            """



### bam files for selected genes
rule run_bam_selected_genes:
    input: "results/mapped/{sample_name}_Aligned.sortedByCoord.out.bam", "results/mapped/{sample_name}_Aligned.sortedByCoord.out.bam.bai"
    output: "results/mapped/{sample_name}_selected_genes.bam", "results/mapped/{sample_name}_selected_genes.bam.bai"
    log:    "00log/run_run_bam_selected_genes_{sample_name}.log"
    conda: "../envs/bioinf_tools.yaml"
    resources: 
        cpu = 2,
        mem = "10",
        time = "34:00:00"
    params: 
        selected_genes = config["selected_genes"]
    message: "run_bam_selected_genes {input}: {resources.cpu} threads / {resources.mem}"
    shell:
        """
        samtools view -b -q 20 -f 3 -L {params.selected_genes}  {input[0]} > {output[0]}
        samtools index {output[0]}
        """
