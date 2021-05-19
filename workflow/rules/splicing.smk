rule majiq_build:
    input: "results/mapped/{sample_name}_Aligned.sortedByCoord.out.bam", "results/mapped/{sample_name}_Aligned.sortedByCoord.out.bam.bai"
    output: "results/splicing/majiq/majiq_{sample_name}/build_{sample_name}/{sample_name}.sj"
    log:  "00log/{sample_name}.bamCoverage"
    #conda: "../envs/deeptools.yaml"
    resources:
        cpu = 4,
        mem = "10G",
        time = "12:00:00"
    params:
        out_folder = "results/splicing/majiq/majiq_{sample_name}",
        settings_file = "results/splicing/majiq/majiq_{sample_name}/settings.txt",
        bamdirs = os.getcwd() + "/results/mapped/",
        sample_name = "{sample_name}",
        gff3 = config['gff3']
    message: "majiq_build {input}: {resources.cpu} threads" #"/ {params.mem}"
    shell: 
        """
        conda activate majiq_env
        [ -d {params.out_folder} ] || mkdir {params.out_folder}
        echo "[info]" >> {params.settings_file}
        echo "readlen=152" >> {params.settings_file}
        echo "bamdirs={params.bamdirs}" >> {params.settings_file}
        echo "genome=hg38" >> {params.settings_file}
        echo "strandness=None" >> {params.settings_file}
        echo "[experiments]" >> {params.settings_file}
        echo "{params.sample_name}={params.sample_name}" >> {params.settings_file}
        cd {params.out_folder}
        majiq build {params.gff3} --nproc {resources.cpu} --junc-files-only -o ./build_{params.sample_name} --simplify --conf {params.settings_file}
        """


rule majiq_build_combine:
    input: "results/mapped/{sample_name}_Aligned.sortedByCoord.out.bam", "results/mapped/{sample_name}_Aligned.sortedByCoord.out.bam.bai"
    output: "results/splicing/majiq/majiq_{sample_name}/build_{sample_name}/{sample_name}.sj"
    log:  "00log/{sample_name}.bamCoverage"
    #conda: "../envs/deeptools.yaml"
    resources:
        cpu = 4,
        mem = "10G",
        time = "12:00:00"
    params:
        out_folder = "results/splicing/majiq/majiq_{sample_name}",
        settings_file = "results/splicing/majiq/majiq_{sample_name}/settings.txt",
        bamdirs = os.getcwd() + "/results/mapped/",
        sample_name = "{sample_name}",
        gff3 = config['gff3']
    message: "majiq_build_combine {input}: {resources.cpu} threads" #"/ {params.mem}"
    shell: 
        """
        conda activate majiq_env
#[experiments]
#SRR1791098=SRR1791098
#SRR1791099=SRR1791099
#SRR1791100=SRR1791100
majiq build /mnt/isilon/thomas-tikhonenko_lab/user/radensc/genomic_info/human/hg38_GRCh38_Reference_Genome/transcriptome_annotation/Homo_sapiens.GRCh38.94.chr.gff3 --nproc 4 --incremental -o ./build --conf settings.txt
        [ -d {params.out_folder} ] || mkdir {params.out_folder}
        echo "[info]" >> {params.settings_file}
        echo "readlen=152" >> {params.settings_file}
        echo "sjdirs={params.bamdirs}" >> {params.settings_file}
        echo "genome=hg38" >> {params.settings_file}
        echo "strandness=None" >> {params.settings_file}
        echo "[experiments]" >> {params.settings_file}
        echo "{params.sample_name}={params.sample_name}" >> {params.settings_file}
        cd {params.out_folder}
        majiq build {params.gff3} --nproc {resources.cpu} --junc-files-only -o ./build_{params.sample_name} --simplify --conf {params.settings_file}
        """
