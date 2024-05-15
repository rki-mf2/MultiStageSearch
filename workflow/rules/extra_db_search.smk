rule ExtraSearchSpectraAgainstProteomes:
    input: 
        proteomes_decoy_fasta = RESULT_DIR / "{sample}/Database/extra_search_proteomes_unique_concatenated_target_decoy.fasta",
        mgf = SearchDB.get_input_MGF()["mgf"],
        par = SEARCHGUI_PAR_FILE,
    output:  
        out_zip = RESULT_DIR / "{sample}/ExtraSearch/proteomes_searchgui_out.zip"
    log:
        stderr_log=RESULT_DIR / "logs/{sample}/ExtraSearch/SearchSpectraAgainstProteomes/stderr.log",
        stdout_log=RESULT_DIR / "logs/{sample}/ExtraSearch/SearchSpectraAgainstProteomes/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/extra_db_search/ExtraSearchSpectraAgainstProteomes.benchmark.txt"
    params:
        searchgui = config["SearchGUI"],
        result_dir = str(RESULT_DIR / "{sample}/ExtraSearch"),
        refname = "proteomes",
        search_engine = "-xtandem",
        psm_fdr = SearchDB.getSearchguiFDRsDBSearch()["psm_fdr"],
        peptide_fdr = SearchDB.getSearchguiFDRsDBSearch()["peptide_fdr"],
        protein_fdr = SearchDB.getSearchguiFDRsDBSearch()["protein_fdr"]
    conda:
        "../envs/java.yml"
    threads: workflow.cores
    resources:
        mem_mb = config["SearchGUI_mem_mb"]
    shell: 
        "java -cp {params.searchgui} eu.isas.searchgui.cmd.SearchCLI -spectrum_files {input.mgf} -fasta_file {input.proteomes_decoy_fasta} -output_folder {params.result_dir} -id_params {input.par} -output_default_name {params.refname}_searchgui_out -psm_fdr {params.psm_fdr} -peptide_fdr {params.peptide_fdr} -protein_fdr {params.protein_fdr} {params.search_engine} 1 -threads {threads} > {log.stdout_log} 2> {log.stderr_log}"


rule extractSearchguiZipExtraSearch:
    input:
        searchgui_zip = RESULT_DIR / "{sample}/ExtraSearch/proteomes_searchgui_out.zip"
    output:
        searchgui_unpacked = directory(RESULT_DIR / "{sample}/ExtraSearch/searchgui_out")
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/ExtraSearch/extractSearchguiZipExtraSearch/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/ExtraSearch/extractSearchguiZipExtraSearch/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/ExtraSearch/extractSearchguiZipExtraSearch.benchmark.txt"
    threads: 1
    shell:
        "unzip {input.searchgui_zip} -d {output.searchgui_unpacked} > {log.stdout_log}"


rule extractXTandemXMLExtraSearch:
    input:
        searchgui_unpacked = RESULT_DIR / "{sample}/ExtraSearch/searchgui_out"
    output:
        XTandem_xml =  RESULT_DIR / "{sample}/ExtraSearch/{sample}.t.xml"
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/ExtraSearch/extractXTandemXMLExtraSearch/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/ExtraSearch/extractXTandemXMLExtraSearch/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/ExtraSearch/extractXTandemXMLExtraSearch.benchmark.txt"
    params:
        XTandem_gz =  str(RESULT_DIR / "{sample}/ExtraSearch/searchgui_out/{sample}.t.xml.gz")
    threads: 1
    shell:
        "gunzip -c {params.XTandem_gz} > {output.XTandem_xml} 2> {log.stderr_log}"


rule runMS2RescoreExtraSearch:
    input:
        XTandem_xml =  RESULT_DIR / "{sample}/ExtraSearch/{sample}.t.xml",
        ms2_config = RESULT_DIR / "{sample}/MS2Rescore/ms2_config.json",
        mgf = MGF_FILE,
        proteomes_decoy_fasta = RESULT_DIR / "{sample}/Database/extra_search_proteomes_unique_concatenated_target_decoy.fasta",
    output: 
        ms2_tsv = RESULT_DIR / "{sample}/ExtraSearch/MS2Rescore/rescored.psms.tsv"
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/ExtraSearch/runMS2RescoreExtraSearch/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/ExtraSearch/runMS2RescoreExtraSearch/stdout.log",
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/ExtraSearch/runMS2RescoreExtraSearch.benchmark.txt"
    params:
        output_path = str(RESULT_DIR / "{sample}/ExtraSearch/MS2Rescore/rescored"),
    conda:
        "../envs/MS2Rescore.yml"
    threads: workflow.cores
    shell:
        "ms2rescore -c {input.ms2_config} -s {input.mgf} -o {params.output_path} -p {input.XTandem_xml} -f {input.proteomes_decoy_fasta} > {log.stdout_log} 2> {log.stderr_log}"


rule convertMS2RescoreTSVExtraSearch:
    input:
        ms2_tsv = RESULT_DIR / "{sample}/ExtraSearch/MS2Rescore/rescored.psms.tsv",
        XTandem_xml =  RESULT_DIR / "{sample}/ExtraSearch/{sample}.t.xml",
        orf_list_mapping = RESULT_DIR / "{sample}/ExtraSearch/ofr_list_id_mapping.tsv",
    output:
        converted_ms2_tsv = RESULT_DIR / "{sample}/ExtraSearch/converted_ms2rescored.tsv"
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/ExtraSearch/convertMS2RescoreTSVExtraSearch/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/ExtraSearch/convertMS2RescoreTSVExtraSearch/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/ExtraSearch/convertMS2RescoreTSVExtraSearch.benchmark.txt"
    params:
        fdr = config["db_search"]["MS2_params"]["fdr"]
    conda:
        "../envs/pandas.yml"
    threads: 1
    script:
        "../scripts/convertMS2rescoredTsv.py"


rule ExtraRunPeptideShakerProteomes:
    input:
        searchgui_zip = RESULT_DIR / "{sample}/ExtraSearch/proteomes_searchgui_out.zip",
        mgf = SearchDB.get_input_MGF()["mgf"],
        proteomes_decoy_fasta = RESULT_DIR / "{sample}/Database/extra_search_proteomes_unique_concatenated_target_decoy.fasta",
    output: 
        peptide_shaker_psdb = RESULT_DIR / "{sample}/ExtraSearch/proteomes.psdb"
    log:
        stderr_log=RESULT_DIR / "logs/{sample}/ExtraSearch/RunPeptideShakerProteomes/stderr.log",
        stdout_log=RESULT_DIR / "logs/{sample}/ExtraSearch/RunPeptideShakerProteomes/stdout.log",
        peptide_shaker_log=RESULT_DIR / "logs/{sample}/ExtraSearch/RunPeptideShakerProteomes/PeptideShaker.log",
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/extra_db_search/ExtraRunPeptideShakerProteomes.benchmark.txt"
    params:
        refname = "proteomes",
        peptideshaker = config["PeptideShaker"],
    conda:
        "../envs/java.yml"
    threads: workflow.cores
    resources:
        mem_mb = config["PeptideShaker_mem_mb"]
    retries: 3 # sometimes there are java exceptions
    shell:
        "java -cp {params.peptideshaker} eu.isas.peptideshaker.cmd.PeptideShakerCLI -reference {params.refname} -fasta_file {input.proteomes_decoy_fasta} -identification_files {input.searchgui_zip} -spectrum_files {input.mgf} -out {output.peptide_shaker_psdb} -threads {threads} -log {log.peptide_shaker_log} > {log.stdout_log} 2> {log.stderr_log}" 


rule ExtraSimplePeptideListProteomes:
    input:
        peptide_shaker_psdb = RESULT_DIR / "{sample}/ExtraSearch/proteomes.psdb"
    output: 
        peptide_shaker_report = RESULT_DIR / "{sample}/ExtraSearch/proteomes_Default_PSM_Report.txt",
    log:
        stderr_log=RESULT_DIR / "logs/{sample}/ExtraSearch/SimplePeptideListProteomes/stderr.log",
        stdout_log=RESULT_DIR / "logs/{sample}/ExtraSearch/SimplePeptideListProteomes/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/extra_db_search/ExtraSimplePeptideListProteomes.benchmark.txt"
    params:
        out_dir = str(RESULT_DIR / "{sample}/ExtraSearch"),
        peptideshaker = config["PeptideShaker"],
    conda:
        "../envs/java.yml"
    threads: workflow.cores
    resources:
        mem_mb = config["PeptideShaker_mem_mb"]
    shell:
        "java -cp {params.peptideshaker} eu.isas.peptideshaker.cmd.ReportCLI -in {input.peptide_shaker_psdb} -out_reports {params.out_dir} -reports 3 > {log.stdout_log} 2> {log.stderr_log}" 


rule ExtraMapORFListsProteomes:
    input:
        peptide_shaker_report = RESULT_DIR / "{sample}/ExtraSearch/proteomes_Default_PSM_Report.txt",
        orf_list_mapping = RESULT_DIR / "{sample}/ExtraSearch/ofr_list_id_mapping.tsv",
    output: 
        converted_peptide_shaker_report = RESULT_DIR / "{sample}/ExtraSearch/converted_proteomes_Default_PSM_Report.txt",
    log:
        stderr_log=RESULT_DIR / "logs/{sample}/ExtraSearch/MapORFListsProteomes/stderr.log",
        stdout_log=RESULT_DIR / "logs/{sample}/ExtraSearch/MapORFListsProteomes/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/extra_db_search/MapORFListsProteomes.benchmark.txt"
    wildcard_constraints:
        sample="[^/]+"
    conda:
        "../envs/pandas.yml"
    threads: workflow.cores
    resources:
        mem_mb = config["PeptideShaker_mem_mb"]
    script:
        "../scripts/convert_peptideshaker_report.py"