rule genome2Proteome:
    input: 
        concat_strain_genomes = CovidMode.checkForCovidMode
    output:
        touch(RESULT_DIR / "{sample}/sixpack/genome2Proteome.done"),
        sixpack_out = RESULT_DIR / "{sample}/sixpack/sixpack_out.txt"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/genome2proteome/genome2Proteome.benchmark.txt"
    params:
        sample_path = str(RESULT_DIR / "{sample}"),
        orfminsize = config["sixpack"]["orfminsize"],
        additional_params = config["sixpack"]["additional_parameters"],
        sixpack_dir_name = "sixpack",
        sixpack_temp = str(RESULT_DIR / "{sample}/sixpack/sixpack_temp.txt"),
        system_latency = config["system_latency"]
    wildcard_constraints:
        sample="[^/]+"
    conda:
        "../envs/sixpack.yml"
    threads: workflow.cores # to ensure no confusion with other samples (happens if the rule runs multiple times at once)
    script:
        "../scripts/genome2proteome.py"

rule concatProteomes:
    input:
        sixpack_done = RESULT_DIR / "{sample}/sixpack/genome2Proteome.done"
    output:
        concat_proteomes = RESULT_DIR / "{sample}/Database/proteomes.fasta"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/genome2proteome/concatProteomes.benchmark.txt"
    params:
        res_dir = RESULT_DIR,
        sample_name = "{sample}",
        sixpack_dir_name = "sixpack"
    wildcard_constraints:
        sample="[^/]+"
    conda:
        "../envs/base_python.yml"
    threads: workflow.cores # to ensure no confusion with other samples (happens if the rule runs multiple times at once)
    script:
        "../scripts/concat_proteomes.py"

rule FilterDuplicateProteomes:
    input:
        strain_accessions = RESULT_DIR / "{sample}/FetchData/strain_accessions.csv",
        concat_proteomes = RESULT_DIR / "{sample}/Database/proteomes.fasta",
    output:
        filtered_proteomes = RESULT_DIR / "{sample}/Database/filtered_proteomes.fasta"
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/genome2proteome/FilterDuplicateProteomes/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/genome2proteome/FilterDuplicateProteomes/stdout.log" 
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/genome2proteome/FilterDuplicateProteomes.benchmark.txt"
    params:
        similarity_threshold = config["similarity_threshold"],
        words_blacklist = config["words_blacklist"],
    wildcard_constraints:
        sample="[^/]+"
    conda:
        "../envs/mapping.yml"
    threads: workflow.cores # to ensure no confusion with other samples (happens if the rule runs multiple times at once)
    script:
        "../scripts/filterDuplicateProteomes.py"

rule AddDecoysProteome:
    input: 
        proteomes = RESULT_DIR / "{sample}/Database/filtered_proteomes.fasta"
    output: 
        proteomes_decoy = RESULT_DIR / "{sample}/Database/filtered_proteomes_concatenated_target_decoy.fasta"
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/genome2proteome/proteomeDB/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/genome2proteome/proteomeDB/stdout.log" 
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/genome2proteome/AddDecoysProteome.benchmark.txt"
    params:
        searchgui = config["SearchGUI"],
    wildcard_constraints:
        sample="[^/]+"
    conda:
        "../envs/java.yml"
    threads: workflow.cores # to ensure no confusion with other samples (happens if the rule runs multiple times at once)
    shell: 
        "java -cp {params.searchgui} eu.isas.searchgui.cmd.FastaCLI -in {input.proteomes} -decoy > {log.stdout_log} 2> {log.stderr_log}"

rule FilterDuplicateORFs:
    input:
        proteomes = RESULT_DIR / "{sample}/Database/filtered_proteomes_concatenated_target_decoy.fasta"
    output:
        proteomes_decoy_unique = RESULT_DIR / "{sample}/Database/proteomes_unique_concatenated_target_decoy.fasta",
        orf_list_mapping = RESULT_DIR / "{sample}/FinalSearch/ofr_list_id_mapping.tsv",
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/genome2proteome/FilterDuplicates/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/genome2proteome/FilterDuplicates/stdout.log" 
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/genome2proteome/FilterDuplicateORFs.benchmark.txt"
    conda:
        "../envs/fetch_data.yml"
    threads: workflow.cores # to ensure no confusion with other samples (happens if the rule runs multiple times at once)
    script:
        "../scripts/filter_and_map_duplicate_ORFs.py"