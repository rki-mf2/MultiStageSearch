rule Extragenome2Proteome:
    input: 
        used_strain_genomes = RESULT_DIR / "{sample}/Phylogeny/used_strain_genomes.fasta",
    output:
        touch(RESULT_DIR / "{sample}/extra_sixpack/genome2Proteome.done"),
        sixpack_out = RESULT_DIR / "{sample}/extra_sixpack/sixpack_out.txt"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/extra_genome2proteome/Extragenome2Proteome.benchmark.txt"
    params:
        sample_path = str(RESULT_DIR / "{sample}"),
        orfminsize = config["sixpack"]["orfminsize"],
        additional_params = config["sixpack"]["additional_parameters"],
        sixpack_dir_name = "extra_sixpack",
        sixpack_temp = str(RESULT_DIR / "{sample}/extra_sixpack/sixpack_temp.txt"),
        system_latency = config["system_latency"]
    wildcard_constraints:
        sample="[^/]+"
    conda:
        "../envs/sixpack.yml"
    threads: workflow.cores # to ensure no confusion with other samples (happens if the rule runs multiple times at once)
    script:
        "../scripts/genome2proteome.py"

rule ExtraconcatProteomes:
    input:
        sixpack_done = RESULT_DIR / "{sample}/extra_sixpack/genome2Proteome.done"
    output:
        concat_proteomes = RESULT_DIR / "{sample}/Database/extra_search_proteomes.fasta"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/extra_genome2proteome/ExtraconcatProteomes.benchmark.txt"
    params:
        res_dir = RESULT_DIR,
        sample_name = "{sample}",
        sixpack_dir_name = "extra_sixpack"
    wildcard_constraints:
        sample="[^/]+"
    conda:
        "../envs/base_python.yml"
    threads: workflow.cores # to ensure no confusion with other samples (happens if the rule runs multiple times at once)
    script:
        "../scripts/concat_proteomes.py"

rule ExtraAddDecoysProteome:
    input: 
        proteomes = RESULT_DIR / "{sample}/Database/extra_search_proteomes.fasta"
    output: 
        proteomes_decoy = RESULT_DIR / "{sample}/Database/extra_search_proteomes_concatenated_target_decoy.fasta"
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/extragenome2proteome/proteomeDB/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/extragenome2proteome/proteomeDB/stdout.log" 
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/extra_genome2proteome/ExtraAddDecoysProteome.benchmark.txt"
    params:
        searchgui = config["SearchGUI"],
    wildcard_constraints:
        sample="[^/]+"
    conda:
        "../envs/java.yml"
    threads: workflow.cores # to ensure no confusion with other samples (happens if the rule runs multiple times at once)
    shell: 
        "java -cp {params.searchgui} eu.isas.searchgui.cmd.FastaCLI -in {input.proteomes} -decoy > {log.stdout_log} 2> {log.stderr_log}"

rule ExtraFilterDuplicateORFs:
    input:
        proteomes = RESULT_DIR / "{sample}/Database/extra_search_proteomes_concatenated_target_decoy.fasta"
    output:
        proteomes_decoy_unique = RESULT_DIR / "{sample}/Database/extra_search_proteomes_unique_concatenated_target_decoy.fasta",
        orf_list_mapping = RESULT_DIR / "{sample}/ExtraSearch/ofr_list_id_mapping.tsv",
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/extragenome2proteome/FilterDuplicates/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/extragenome2proteome/FilterDuplicates/stdout.log" 
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/extra_genome2proteome/ExtraFilterDuplicateORFs.benchmark.txt"
    wildcard_constraints:
        sample="[^/]+"
    conda:
        "../envs/fetch_data.yml"
    threads: 1
    script:
        "../scripts/filter_and_map_duplicate_ORFs.py"