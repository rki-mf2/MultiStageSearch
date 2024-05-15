checkpoint splitLineages:
    input:
        accessions_df = RESULT_DIR / "{sample}/FetchData/strain_accessions.csv",
        concat_fasta = RESULT_DIR / "{sample}/FetchData/concat_strain_genomes.fasta"
    output:
        splitted_fasta_files = directory(RESULT_DIR / "{sample}/CovidMode/Lineages"),
        splitLineages_done = touch(RESULT_DIR / "{sample}/CovidMode/Lineages/splitLineages.done")
    log:
        stderr_log=RESULT_DIR / "logs/{sample}/CovidMode/splitLineages/stderr.log",
        stdout_log=RESULT_DIR / "logs/{sample}/CovidMode/splitLineages/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/CovidMode/splitLineages.benchmark.txt"
    params:
        dir_path = str(RESULT_DIR / "{sample}/CovidMode/Lineages"),
        num_lineages = config["CovidMode"]["min_num_genomes_per_lineage"],
    wildcard_constraints:
        sample="[^/]+"
    conda:
        "../envs/fetch_data.yml"
    threads: workflow.cores
    script:
        "../scripts/splitLineages.py"


rule CovidGenome2Proteome:
    input: 
        dynamic_files = CovidMode.get_input_Covidgenome2Proteome
    output:
        touch(RESULT_DIR / "{sample}/CovidMode/{lineage}/sixpack/genome2Proteome.done"),
        sixpack_out = RESULT_DIR / "{sample}/CovidMode/{lineage}/sixpack/sixpack_out.txt"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/CovidMode/{lineage}/sixpack/genome2Proteome.benchmark.txt"
    params:
        sample_path = str(RESULT_DIR / "{sample}/CovidMode/{lineage}"),
        orfminsize = config["sixpack"]["orfminsize"],
        additional_params = config["sixpack"]["additional_parameters"],
        sixpack_dir_name = "sixpack",
        sixpack_temp = str(RESULT_DIR / "{sample}/CovidMode/{lineage}/sixpack/sixpack_temp.txt"),
        system_latency = config["system_latency"]
    wildcard_constraints:
        sample="[^/]+"
    conda:
        "../envs/sixpack.yml"
    threads: workflow.cores
    script:
        "../scripts/genome2proteome.py"


rule CovidConcatProteomes:
    input:
        sixpack_done = RESULT_DIR / "{sample}/CovidMode/{lineage}/sixpack/genome2Proteome.done"
    output:
        concat_proteomes = RESULT_DIR / "{sample}/CovidMode/{lineage}/proteomes.fasta"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/CovidMode/{lineage}/CovidconcatProteomes.benchmark.txt"
    params:
        res_dir = RESULT_DIR,
        sample_name = "{sample}/CovidMode/{lineage}/",
        sixpack_dir_name = "sixpack"
    wildcard_constraints:
        sample="[^/]+"
    conda:
        "../envs/base_python.yml"
    threads: workflow.cores
    script:
        "../scripts/concat_proteomes.py"


rule CovidAddDecoysProteome:
    input: 
        proteomes = RESULT_DIR / "{sample}/CovidMode/{lineage}/proteomes.fasta"
    output: 
        proteomes_decoy = RESULT_DIR / "{sample}/CovidMode/{lineage}/proteomes_concatenated_target_decoy.fasta"
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/CovidMode/{lineage}/CovidAddDecoysProteome/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/CovidMode/{lineage}/CovidAddDecoysProteome/stdout.log" 
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/CovidMode/{lineage}/CovidAddDecoysProteome.benchmark.txt"
    params:
        searchgui = config["SearchGUI"],
    wildcard_constraints:
        sample="[^/]+"
    conda:
        "../envs/java.yml"
    threads: workflow.cores
    shell: 
        "java -cp {params.searchgui} eu.isas.searchgui.cmd.FastaCLI -in {input.proteomes} -decoy > {log.stdout_log} 2> {log.stderr_log}"


rule CovidFilterDuplicateORFs:
    input:
        proteomes = RESULT_DIR / "{sample}/CovidMode/{lineage}/proteomes_concatenated_target_decoy.fasta"
    output:
        proteomes_decoy_unique = RESULT_DIR / "{sample}/CovidMode/{lineage}/proteomes_unique_concatenated_target_decoy.fasta",
        orf_list_mapping = RESULT_DIR / "{sample}/CovidMode/{lineage}/LineageSearch/ofr_list_id_mapping.tsv",
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/CovidMode/{lineage}/CovidFilterDuplicateORFs/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/CovidMode/{lineage}/CovidFilterDuplicateORFs/stdout.log" 
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/CovidMode/{lineage}/CovidFilterDuplicateORFs.benchmark.txt"
    conda:
        "../envs/fetch_data.yml"
    threads: workflow.cores
    script:
        "../scripts/filter_and_map_duplicate_ORFs.py"


rule CovidSearchSpectraAgainstProteomes:
    input: 
        proteomes_decoy_fasta = RESULT_DIR / "{sample}/CovidMode/{lineage}/proteomes_unique_concatenated_target_decoy.fasta",
        mgf = SearchDB.get_input_MGF()["mgf"],
        par = SEARCHGUI_PAR_FILE,
    output:  
        out_zip = RESULT_DIR / "{sample}/CovidMode/{lineage}/proteomes_searchgui_out.zip"
    log:
        stderr_log=RESULT_DIR / "logs/{sample}/CovidMode/{lineage}/CovidSearchSpectraAgainstProteomes/stderr.log",
        stdout_log=RESULT_DIR / "logs/{sample}/CovidMode/{lineage}/CovidSearchSpectraAgainstProteomes/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/final_db_search/{lineage}/CovidSearchSpectraAgainstProteomes.benchmark.txt"
    params:
        searchgui = config["SearchGUI"],
        result_dir = str(RESULT_DIR / "{sample}/CovidMode/{lineage}"),
        refname = "proteomes",
        search_engine = "-xtandem",
        psm_fdr = SearchDB.getSearchguiFDRsDBSearch()["psm_fdr"],
        peptide_fdr = SearchDB.getSearchguiFDRsDBSearch()["peptide_fdr"],
        protein_fdr = SearchDB.getSearchguiFDRsDBSearch()["protein_fdr"]
    wildcard_constraints:
        sample="[^/]+"
    conda:
        "../envs/java.yml"
    threads: workflow.cores
    resources:
        mem_mb = config["SearchGUI_mem_mb"]
    shell: 
        "java -cp {params.searchgui} eu.isas.searchgui.cmd.SearchCLI -spectrum_files {input.mgf} -fasta_file {input.proteomes_decoy_fasta} -output_folder {params.result_dir} -id_params {input.par} -output_default_name {params.refname}_searchgui_out -psm_fdr {params.psm_fdr} -peptide_fdr {params.peptide_fdr} -protein_fdr {params.protein_fdr} {params.search_engine} 1 -threads {threads} > {log.stdout_log} 2> {log.stderr_log}"


rule CovidExtractSearchguiZipFinalSearch:
    input:
        searchgui_zip = RESULT_DIR / "{sample}/CovidMode/{lineage}/proteomes_searchgui_out.zip"
    output:
        searchgui_unpacked = directory(RESULT_DIR / "{sample}/CovidMode/{lineage}/searchgui_out")
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/CovidMode/{lineage}/CovidExtractSearchguiZipFinalSearch/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/CovidMode/{lineage}/CovidExtractSearchguiZipFinalSearch/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/CovidMode/{lineage}/CovidExtractSearchguiZipFinalSearch.benchmark.txt"
    wildcard_constraints:
        sample="[^/]+"
    threads: 1
    shell:
        "unzip {input.searchgui_zip} -d {output.searchgui_unpacked} > {log.stdout_log}"


rule CovidExtractXTandemXMLFinalSearch:
    input:
        searchgui_unpacked = RESULT_DIR / "{sample}/CovidMode/{lineage}/searchgui_out"
    output:
        XTandem_xml =  RESULT_DIR / "{sample}/CovidMode/{lineage}/{sample}.t.xml"
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/CovidMode/{lineage}/CovidExtractXTandemXMLFinalSearch/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/CovidMode/{lineage}/CovidExtractXTandemXMLFinalSearch/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/CovidMode/{lineage}/CovidExtractXTandemXMLFinalSearch.benchmark.txt"
    params:
        XTandem_gz =  str(RESULT_DIR / "{sample}/CovidMode/{lineage}/searchgui_out/{sample}.t.xml.gz")
    wildcard_constraints:
        sample="[^/]+"
    threads: 1
    shell:
        "gunzip -c {params.XTandem_gz} > {output.XTandem_xml} 2> {log.stderr_log}"


rule CovidRunMS2RescoreFinalSearch:
    input:
        XTandem_xml =  RESULT_DIR / "{sample}/CovidMode/{lineage}/{sample}.t.xml",
        ms2_config = RESULT_DIR / "{sample}/MS2Rescore/ms2_config.json",
        mgf = MGF_FILE,
        proteomes_decoy_fasta = RESULT_DIR / "{sample}/CovidMode/{lineage}/proteomes_unique_concatenated_target_decoy.fasta",
    output: 
        ms2_tsv = RESULT_DIR / "{sample}/CovidMode/{lineage}/MS2Rescore/rescored.psms.tsv"
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/CovidMode/{lineage}/CovidRunMS2RescoreFinalSearch/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/CovidMode/{lineage}/CovidRunMS2RescoreFinalSearch/stdout.log",
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/CovidMode/{lineage}/CovidRunMS2RescoreFinalSearch.benchmark.txt"
    params:
        output_path = str(RESULT_DIR / "{sample}/CovidMode/{lineage}/MS2Rescore/rescored"),
    wildcard_constraints:
        sample="[^/]+"
    conda:
        "../envs/MS2Rescore.yml"
    threads: workflow.cores
    shell:
        "ms2rescore -c {input.ms2_config} -s {input.mgf} -o {params.output_path} -p {input.XTandem_xml} -f {input.proteomes_decoy_fasta} > {log.stdout_log} 2> {log.stderr_log}"


rule CovidConvertMS2RescoreTSVFinalSearch:
    input:
        ms2_tsv = RESULT_DIR / "{sample}/CovidMode/{lineage}/MS2Rescore/rescored.psms.tsv",
        XTandem_xml =  RESULT_DIR / "{sample}/CovidMode/{lineage}/{sample}.t.xml",
        orf_list_mapping = RESULT_DIR / "{sample}/CovidMode/{lineage}/LineageSearch/ofr_list_id_mapping.tsv",
    output:
        converted_ms2_tsv = RESULT_DIR / "{sample}/CovidMode/{lineage}/converted_ms2rescored.tsv"
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/CovidMode/{lineage}/CovidConvertMS2RescoreTSVFinalSearch/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/CovidMode/{lineage}/CovidConvertMS2RescoreTSVFinalSearch/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/CovidMode/{lineage}/CovidConvertMS2RescoreTSVFinalSearch.benchmark.txt"
    params:
        fdr = config["db_search"]["MS2_params"]["fdr"]
    wildcard_constraints:
        sample="[^/]+"
    conda:
        "../envs/pandas.yml"
    threads: 1
    script:
        "../scripts/convertMS2rescoredTsv.py"


rule CovidRunPeptideShakerProteomes:
    input:
        searchgui_zip = RESULT_DIR / "{sample}/CovidMode/{lineage}/proteomes_searchgui_out.zip",
        mgf = SearchDB.get_input_MGF()["mgf"],
        proteomes_decoy_fasta = RESULT_DIR / "{sample}/CovidMode/{lineage}/proteomes_unique_concatenated_target_decoy.fasta",
    output: 
        peptide_shaker_psdb = RESULT_DIR / "{sample}/CovidMode/{lineage}/proteomes.psdb"
    log:
        stderr_log=RESULT_DIR / "logs/{sample}/CovidMode/{lineage}/CovidRunPeptideShakerProteomes/stderr.log",
        stdout_log=RESULT_DIR / "logs/{sample}/CovidMode/{lineage}/CovidRunPeptideShakerProteomes/stdout.log",
        peptide_shaker_log=RESULT_DIR / "logs/{sample}/CovidMode/{lineage}/CovidRunPeptideShakerProteomes/PeptideShaker.log",
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/CovidMode/{lineage}/CovidRunPeptideShakerProteomes.benchmark.txt"
    params:
        refname = "proteomes",
        peptideshaker = config["PeptideShaker"],
    wildcard_constraints:
        sample="[^/]+"
    conda:
        "../envs/java.yml"
    threads: workflow.cores
    resources:
        mem_mb = config["PeptideShaker_mem_mb"]
    retries: 3 # sometimes there are java exceptions
    shell:
        "java -cp {params.peptideshaker} eu.isas.peptideshaker.cmd.PeptideShakerCLI -reference {params.refname} -fasta_file {input.proteomes_decoy_fasta} -identification_files {input.searchgui_zip} -spectrum_files {input.mgf} -out {output.peptide_shaker_psdb} -threads {threads} -log {log.peptide_shaker_log} > {log.stdout_log} 2> {log.stderr_log}" 


rule CovidSimplePeptideListProteomes:
    input:
        peptide_shaker_psdb = RESULT_DIR / "{sample}/CovidMode/{lineage}/proteomes.psdb"
    output: 
        peptide_shaker_report = RESULT_DIR / "{sample}/CovidMode/{lineage}/proteomes_Default_PSM_Report.txt",
    log:
        stderr_log=RESULT_DIR / "logs/{sample}/CovidMode/{lineage}/SimplePeptideListProteomes/stderr.log",
        stdout_log=RESULT_DIR / "logs/{sample}/CovidMode/{lineage}/SimplePeptideListProteomes/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/CovidMode/{lineage}/SimplePeptideListProteomes.benchmark.txt"
    params:
        out_dir = str(RESULT_DIR / "{sample}/CovidMode/{lineage}"),
        peptideshaker = config["PeptideShaker"],
    wildcard_constraints:
        sample="[^/]+"
    conda:
        "../envs/java.yml"
    threads: workflow.cores
    resources:
        mem_mb = config["PeptideShaker_mem_mb"]
    shell:
        "java -cp {params.peptideshaker} eu.isas.peptideshaker.cmd.ReportCLI -in {input.peptide_shaker_psdb} -out_reports {params.out_dir} -reports 3 > {log.stdout_log} 2> {log.stderr_log}" 


rule CovidMapORFListsProteomes:
    input:
        peptide_shaker_report = RESULT_DIR / "{sample}/CovidMode/{lineage}/proteomes_Default_PSM_Report.txt",
        orf_list_mapping = RESULT_DIR / "{sample}/CovidMode/{lineage}/LineageSearch/ofr_list_id_mapping.tsv",
    output: 
        converted_peptide_shaker_report = RESULT_DIR / "{sample}/CovidMode/{lineage}/converted_proteomes_Default_PSM_Report.txt",
    log:
        stderr_log=RESULT_DIR / "logs/{sample}/CovidMode/{lineage}/CovidMapORFListsProteomes/stderr.log",
        stdout_log=RESULT_DIR / "logs/{sample}/CovidMode/{lineage}/CovidMapORFListsProteomes/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/CovidMode/{lineage}/CovidMapORFListsProteomes.benchmark.txt"
    wildcard_constraints:
        sample="[^/]+"
    conda:
        "../envs/pandas.yml"
    threads: workflow.cores
    resources:
        mem_mb = config["PeptideShaker_mem_mb"]
    script:
        "../scripts/convert_peptideshaker_report.py"


rule CovidGetLineageNames:
    input:
        report = SearchDB.getReports()["covid_search_report"],
        strain_accessions = RESULT_DIR / "{sample}/FetchData/strain_accessions.csv",
        fasta = RESULT_DIR / "{sample}/CovidMode/{lineage}/proteomes.fasta"
    output:
        lineage_name_counts = RESULT_DIR / "{sample}/CovidMode/{lineage}/taxids/lineage_name_counts.tsv"
    log:
        stderr_log=RESULT_DIR / "logs/{sample}/CovidMode/{lineage}/FetchData/CovidGetLineageNames/stderr.log",
        stdout_log=RESULT_DIR / "logs/{sample}/CovidMode/{lineage}/FetchData/CovidGetLineageNames/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/CovidMode/{lineage}/fetch_data/CovidGetLineageNames.benchmark.txt"
    params:
        number_taxids = config["mapping"]["number_of_strains"],
    wildcard_constraints:
        sample="[^/]+"
    conda:
        "../envs/fetch_data.yml"
    threads: 1
    script:
        "../scripts/get_species_strain.py"



rule CovidMapORFIDs:
    input:
        fasta = RESULT_DIR / "{sample}/CovidMode/{lineage}/proteomes_concatenated_target_decoy.fasta",
        report = SearchDB.getReports()["covid_search_report"],
        lineage_name_counts = RESULT_DIR / "{sample}/CovidMode/{lineage}/taxids/lineage_name_counts.tsv",
    output:
        all_mapped_ORFs = RESULT_DIR / "{sample}/CovidMode/{lineage}/taxids/mapped_ORFs.txt",
        ORF_scores = RESULT_DIR / "{sample}/CovidMode/{lineage}/taxids/ORF_scores.csv",
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/CovidMode/{lineage}/taxids/CovidGetLineageNames/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/CovidMode/{lineage}/taxids/CovidGetLineageNames/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/CovidMode/{lineage}/mapping/CovidGetLineageNames.benchmark.txt"
    wildcard_constraints:
        sample="[^/]+"
    conda: 
        "../envs/mapping.yml"
    script: 
        "../scripts/ORFMapping.py"


rule CovidAggregateLineageResults:
    input:
        dynamic_files = CovidMode.get_Covid_topScoring,
        concat_fasta = RESULT_DIR / "{sample}/FetchData/concat_strain_genomes.fasta"
    output:
        top_scoring_fasta = RESULT_DIR / "{sample}/CovidMode/top_scoring_lineages.fasta",
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/CovidMode/CovidAggregateLineageResults/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/CovidMode/CovidAggregateLineageResults/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/CovidMode/CovidAggregateLineageResults.benchmark.txt"
    params:
        n_top_scoring = config["CovidMode"]["num_top_scoring_per_lineage"]
    wildcard_constraints:
        sample="[^/]+"
    conda: 
        "../envs/mapping.yml"
    script: 
        "../scripts/aggregateCovidGenomes.py"