rule run_convert_SJ_to_interact:
    input: "results/mapped/{sample_name}_Aligned.sortedByCoord.out.bam", "results/mapped/{sample_name}_Aligned.sortedByCoord.out.bam.bai"
    output: "results/junctions/{sample_name}.interact"
    log:    "00log/run_convert_SJ_to_interact_{sample_name}.log"
    resources: 
        cpu = 2,
        mem = "40",
        time = "44:00:00"
    params: 
        infile = "results/mapped/{sample_name}_SJ.out.tab"
    message: "run_convert_SJ_to_interact {params}: {resources.cpu} threads / {resources.mem}"
    shell:
        "ruby {workflow.basedir}/scripts/convert_SJ_to_interaction.rb {params.infile} {output} 5"

rule run_convert_interact_to_bigbed:
    input: "results/junctions/{sample_name}.interact"
    output: "results/junctions/{sample_name}.inter.bb"
    log:    "00log/run_convert_interact_to_bigbed_{sample_name}.log"
    resources: 
        cpu = 2,
        mem = "40",
        time = "44:00:00"
    params: 
        chrom_sizes = config["chrom_sizes"]
    message: "run_convert_SJ_to_interact {params}: {resources.cpu} threads / {resources.mem}"
    shell:
        """
        grep ^chr {input} > {input}.tmp
        LC_COLLATE=C sort -k1,1 -k2,2n {input}.tmp > {input}.sorted
        # this interact.as file is from ucsc
        bedToBigBed -as={workflow.basedir}/../resources/interact.as  -type=bed5+13 {input}.sorted {params.chrom_sizes} {output}
        rm {input}.tmp
        rm {input}.sorted
        """


rule majiq_build:
    input: "results/mapped/{sample_name}_Aligned.sortedByCoord.out.bam", "results/mapped/{sample_name}_Aligned.sortedByCoord.out.bam.bai"
    output: "results/splicing/majiq/majiq_{sample_name}/build_{sample_name}/{sample_name}.sj"
    log:  "00log/{sample_name}.majiq_build"
    #conda: "../envs/deeptools.yaml"
    resources:
        cpu = 4,
        mem = "10",
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
        [ -d {params.out_folder} ] || mkdir {params.out_folder}
        echo "[info]" > {params.settings_file}
        echo "readlen=152" >> {params.settings_file}
        echo "bamdirs={params.bamdirs}" >> {params.settings_file}
        echo "genome=hg38" >> {params.settings_file}
        echo "strandness=None" >> {params.settings_file}
        echo "[experiments]" >> {params.settings_file}
        echo "{params.sample_name}={params.sample_name}_Aligned.sortedByCoord.out" >> {params.settings_file}
        cd {params.out_folder}
        bash -c '
            . $HOME/.bashrc # if not loaded automatically
            conda activate majiq_env
            majiq build {params.gff3} --nproc {resources.cpu} --junc-files-only -o ./build_{params.sample_name} --simplify --conf settings.txt
            mv build_{params.sample_name}/{params.sample_name}_Aligned.sortedByCoord.out.sj build_{params.sample_name}/{params.sample_name}.sj
            conda deactivate'
        
        """




#rule majiq_build_combine:
#    input: get_sj_files
#    output: "results/splicing/majiq/build"
#    log:  "00log/{sample_name}.majiq_build_combine"
#    #conda: "../envs/deeptools.yaml"
#    resources:
#        cpu = 4,
#        mem = "10G",
#        time = "12:00:00"
#    params:
#        out_folder = "results/splicing/majiq/majiq_{sample_name}",
#        settings_file = "results/splicing/majiq/majiq_{sample_name}/settings.txt",
#        bamdirs = os.getcwd() + "/results/mapped/",
#        sample_name = "{sample_name}",
#        gff3 = config['gff3']
#    message: "majiq_build_combine {input}: {resources.cpu} threads" #"/ {params.mem}"
#    shell: 
#        """
#        [ -d {params.out_folder} ] || mkdir {params.out_folder}
#        echo "[info]" > {params.settings_file}
#        echo "readlen=152" >> {params.settings_file}
#        echo "bamdirs={params.bamdirs}" >> {params.settings_file}
#        echo "genome=hg38" >> {params.settings_file}
#        echo "strandness=None" >> {params.settings_file}
#        echo "[experiments]" >> {params.settings_file}
##[experiments]
##SRR1791098=SRR1791098
##SRR1791099=SRR1791099
##SRR1791100=SRR1791100
#majiq build --incremental -o ./build --conf settings.txt
#        [ -d {params.out_folder} ] || mkdir {params.out_folder}
#        echo "[info]" >> {params.settings_file}
#        echo "readlen=152" >> {params.settings_file}
#        echo "sjdirs={params.bamdirs}" >> {params.settings_file}
#        echo "genome=hg38" >> {params.settings_file}
#        echo "strandness=None" >> {params.settings_file}
#        echo "[experiments]" >> {params.settings_file}
#        echo "{params.sample_name}={params.sample_name}" >> {params.settings_file}
#        cd {params.out_folder}
#        majiq build {params.gff3} --incremental -o ./build --conf settings.txt --nproc {resources.cpu} --simplify
#        """
#