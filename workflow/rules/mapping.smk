rule mapTaxIDs:
    input:
        report = SearchDB.getReports()["first_search_report"],
        taxIDs = config["db_search"]["protein_accessions"],
    output:
        all_mapped_taxIDs = RESULT_DIR / "{sample}/taxids/mapped_taxids.txt",
        taxID_scores = RESULT_DIR / "{sample}/taxids/taxid_scores.tsv",
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/taxids/mapTaxIDs/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/taxids/mapTaxIDs/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/mapping/mapTaxIDs.benchmark.txt"
    wildcard_constraints:
        sample="[^/]+"
    conda: 
        "../envs/mapping.yml"
    script: 
        "../scripts/taxIDMapping.py"


rule mapORFIDs:
    input:
        fasta = RESULT_DIR / "{sample}/Database/filtered_proteomes_concatenated_target_decoy.fasta",
        report = SearchDB.getReports()["final_search_report"],
        strain_taxids = RESULT_DIR / "{sample}/taxids/strain_name_counts.tsv",
    output:
        all_mapped_ORFs = RESULT_DIR / "{sample}/taxids/mapped_ORFs.txt",
        ORF_scores = RESULT_DIR / "{sample}/taxids/ORF_scores.csv",
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/taxids/mapORFIDs/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/taxids/mapORFIDs/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/mapping/mapORFIDs.benchmark.txt"
    wildcard_constraints:
        sample="[^/]+"
    conda: 
        "../envs/mapping.yml"
    script: 
        "../scripts/ORFMapping.py"


rule extraMapORFIDs:
    input:
        fasta = RESULT_DIR / "{sample}/Database/extra_search_proteomes_concatenated_target_decoy.fasta",
        report = SearchDB.getReports()["extra_search_report"],
        strain_taxids = RESULT_DIR / "{sample}/taxids/extra_search_strain_name_counts.tsv",
    output:
        all_mapped_ORFs = RESULT_DIR / "{sample}/taxids/extra_search_mapped_ORFs.txt",
        ORF_scores = RESULT_DIR / "{sample}/taxids/extra_search_ORF_scores.csv",
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/taxids/extra_search_mapORFIDs/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/taxids/extra_search_mapORFIDs/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/mapping/extraMapORFIDs.benchmark.txt"
    wildcard_constraints:
        sample="[^/]+"
    conda: 
        "../envs/mapping.yml"
    script: 
        "../scripts/ORFMapping.py"