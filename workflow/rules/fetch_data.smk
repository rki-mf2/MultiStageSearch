checkpoint fetchStrainGenomes:
    input:
        taxID_scores = RESULT_DIR / "{sample}/taxids/taxid_scores.tsv",
    output:
        out_file_df = RESULT_DIR / "{sample}/FetchData/strain_accessions.csv",
        concat_fasta = RESULT_DIR / "{sample}/FetchData/concat_strain_genomes.fasta",
        check_covid_mode = RESULT_DIR / "{sample}/FetchData/activate_covid_mode.txt"
    log:
        stderr_log=RESULT_DIR / "logs/{sample}/FetchData/fetchStrainGenomes/stderr.log",
        stdout_log=RESULT_DIR / "logs/{sample}/FetchData/fetchStrainGenomes/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/fetch_data/fetchStrainGenomes.benchmark.txt"
    params:
        APIMail = config["Entrez"]["APIMail"],
        APIKey = config["Entrez"]["APIKey"],
        number_taxids = config["mapping"]["number_of_taxids"],
        weight_diff = config["mapping"]["max_weight_differences"],
        max_sequence_length = config["fetchData"]["max_sequence_length"],
        sequence_length_diff = config["fetchData"]["sequence_length_diff"],
        max_number_accessions = config["fetchData"]["max_number_accessions"],
        NCBI = config["fetchData"]["use_NCBI_Taxa"],
        only_use_complete_genomes = config["fetchData"]["only_use_complete_genomes"],
        sqlite_db_path = config["fetchData"]["sqlite_db_path"],
        out_path = str(RESULT_DIR / "{sample}/FetchData"),
    wildcard_constraints:
        sample="[^/]+"
    conda:
        "../envs/fetch_data.yml"
    threads: 1
    retries: 3
    script:
        "../scripts/fetch_data.py"


rule get_strain_names:
    input:
        report = SearchDB.getReports()["final_search_report"],
        strain_accessions = RESULT_DIR / "{sample}/FetchData/strain_accessions.csv",
        proteomes_fasta = RESULT_DIR / "{sample}/Database/filtered_proteomes.fasta",
    output:
        strain_name_counts = RESULT_DIR / "{sample}/taxids/strain_name_counts.tsv"
    log:
        stderr_log=RESULT_DIR / "logs/{sample}/FetchData/get_strain_names/stderr.log",
        stdout_log=RESULT_DIR / "logs/{sample}/FetchData/get_strain_names/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/fetch_data/get_strain_names.benchmark.txt"
    params:
        number_taxids = config["mapping"]["number_of_strains"],
    wildcard_constraints:
        sample="[^/]+"
    conda:
        "../envs/fetch_data.yml"
    threads: 1
    script:
        "../scripts/get_species_strain.py"


rule extra_get_strain_names:
    input:
        report = SearchDB.getReports()["extra_search_report"],
        strain_accessions = RESULT_DIR / "{sample}/FetchData/strain_accessions.csv",
        proteomes_fasta = RESULT_DIR / "{sample}/Database/extra_search_proteomes.fasta",
    output:
        strain_name_counts = RESULT_DIR / "{sample}/taxids/extra_search_strain_name_counts.tsv"
    log:
        stderr_log=RESULT_DIR / "logs/{sample}/FetchData/extra_get_strain_names/stderr.log",
        stdout_log=RESULT_DIR / "logs/{sample}/FetchData/extra_get_strain_names/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/fetch_data/extra_get_strain_names.benchmark.txt"
    params:
        number_taxids = config["mapping"]["number_of_strains"],
    wildcard_constraints:
        sample="[^/]+"
    conda:
        "../envs/fetch_data.yml"
    threads: 1
    script:
        "../scripts/get_species_strain.py"
