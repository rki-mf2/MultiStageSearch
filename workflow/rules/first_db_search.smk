rule AddDecoysRef:
    input: 
        ref = config["db_search"]["ref"],
    output: 
        ref_decoy_fasta = SearchDB.get_output_AddDecoysRef()["ref_decoy_fasta"]
    log:
        stderr_log=RESULT_DIR / "logs/FirstSearch/RefDB/stderr.log",
        stdout_log=RESULT_DIR / "logs/FirstSearch/RefDB/stdout.log" 
    benchmark:
        RESULT_DIR/ "benchmarks/first_db_search/AddDecoysRef.benchmark.txt"
    params:
        searchgui = config["SearchGUI"],
    conda:
        "../envs/java.yml"
    threads: 1
    shell: 
        "java -cp {params.searchgui} eu.isas.searchgui.cmd.FastaCLI -in {input.ref} -decoy > {log.stdout_log} 2> {log.stderr_log}"


rule SearchSpectraAgainstReference:
    input: 
        ref_decoy_fasta = SearchDB.get_output_AddDecoysRef()["ref_decoy_fasta"],
        mgf = SearchDB.get_input_MGF()["mgf"],
        par = SEARCHGUI_PAR_FILE,
    output:  
        out_zip = RESULT_DIR / "{sample}/FirstSearch/ref_searchgui_out.zip"
    log:
        stderr_log=RESULT_DIR / "logs/{sample}/FirstSearch/SearchSpectraAgainstReference/stderr.log",
        stdout_log=RESULT_DIR / "logs/{sample}/FirstSearch/SearchSpectraAgainstReference/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/first_db_search/SearchSpectraAgainstReference.benchmark.txt"
    params:
        searchgui = config["SearchGUI"],
        result_dir = str(RESULT_DIR / "{sample}/FirstSearch"),
        refname = "ref",
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
        "java -cp {params.searchgui} eu.isas.searchgui.cmd.SearchCLI -spectrum_files {input.mgf} -fasta_file {input.ref_decoy_fasta} -output_folder {params.result_dir} -id_params {input.par} -output_default_name {params.refname}_searchgui_out -psm_fdr {params.psm_fdr} -peptide_fdr {params.peptide_fdr} -protein_fdr {params.protein_fdr} {params.search_engine} 1 -threads {threads} > {log.stdout_log} 2> {log.stderr_log}"


rule extractSearchguiZipFirstSearch:
    input:
        searchgui_zip = RESULT_DIR / "{sample}/FirstSearch/ref_searchgui_out.zip"
    output:
        searchgui_unpacked = directory(RESULT_DIR / "{sample}/FirstSearch/searchgui_out")
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/FirstSearch/extractSearchguiZipFirstSearch/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/FirstSearch/extractSearchguiZipFirstSearch/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/first_db_search/extractSearchguiZipFirstSearch.benchmark.txt"
    wildcard_constraints:
        sample="[^/]+"
    threads: 1
    shell:
        "unzip {input.searchgui_zip} -d {output.searchgui_unpacked} > {log.stdout_log}"


rule extractXTandemXMLFirstSearch:
    input:
        searchgui_unpacked = RESULT_DIR / "{sample}/FirstSearch/searchgui_out"
    output:
        XTandem_xml =  RESULT_DIR / "{sample}/FirstSearch/{sample}.t.xml"
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/FirstSearch/extractXTandemXMLFirstSearch/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/FirstSearch/extractXTandemXMLFirstSearch/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/FirstSearch/extractXTandemXMLFirstSearch.benchmark.txt"
    params:
        XTandem_gz =  str(RESULT_DIR / "{sample}/FirstSearch/searchgui_out/{sample}.t.xml.gz")
    wildcard_constraints:
        sample="[^/]+"
    threads: 1
    shell:
        "gunzip -c {params.XTandem_gz} > {output.XTandem_xml} 2> {log.stderr_log}"


rule runMS2RescoreFirstSearch:
    input:
        XTandem_xml =  RESULT_DIR / "{sample}/FirstSearch/{sample}.t.xml",
        ms2_config = RESULT_DIR / "{sample}/MS2Rescore/ms2_config.json",
        mgf = MGF_FILE,
        ref_decoy_fasta = SearchDB.get_output_AddDecoysRef()["ref_decoy_fasta"],
    output: 
        ms2_tsv = RESULT_DIR / "{sample}/FirstSearch/MS2Rescore/rescored.psms.tsv"
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/FirstSearch/runMS2RescoreFirstSearch/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/FirstSearch/runMS2RescoreFirstSearch/stdout.log",
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/FirstSearch/runMS2RescoreFirstSearch.benchmark.txt"
    params:
        output_path = str(RESULT_DIR / "{sample}/FirstSearch/MS2Rescore/rescored"),
    wildcard_constraints:
        sample="[^/]+"
    conda:
        "../envs/MS2Rescore.yml"
    threads: workflow.cores
    shell:
        "ms2rescore -c {input.ms2_config} -s {input.mgf} -o {params.output_path} -p {input.XTandem_xml} -f {input.ref_decoy_fasta} > {log.stdout_log} 2> {log.stderr_log}"


rule convertMS2RescoreTSVFirstSearch:
    input:
        ms2_tsv = RESULT_DIR / "{sample}/FirstSearch/MS2Rescore/rescored.psms.tsv",
        XTandem_xml =  RESULT_DIR / "{sample}/FirstSearch/{sample}.t.xml"
    output:
        converted_ms2_tsv = RESULT_DIR / "{sample}/FirstSearch/converted_ms2rescored.tsv"
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/FirstSearch/convertMS2RescoreTSVFirstSearch/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/FirstSearch/convertMS2RescoreTSVFirstSearch/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/FirstSearch/convertMS2RescoreTSVFirstSearch.benchmark.txt"
    params:
        fdr = config["db_search"]["MS2_params"]["fdr"]
    wildcard_constraints:
        sample="[^/]+"
    conda:
        "../envs/pandas.yml"
    threads: 1
    script:
        "../scripts/convertMS2rescoredTsv.py"


rule RunPeptideShakerRef:
    input:
        searchgui_zip = RESULT_DIR / "{sample}/FirstSearch/ref_searchgui_out.zip",
        mgf = SearchDB.get_input_MGF()["mgf"],
        ref_decoy_fasta = SearchDB.get_output_AddDecoysRef()["ref_decoy_fasta"],
    output: 
        peptide_shaker_psdb = RESULT_DIR / "{sample}/FirstSearch/ref.psdb"
    log:
        stderr_log=RESULT_DIR / "logs/{sample}/FirstSearch/RunPeptideShakerRef/stderr.log",
        stdout_log=RESULT_DIR / "logs/{sample}/FirstSearch/RunPeptideShakerRef/stdout.log",
        peptide_shaker_log=RESULT_DIR / "logs/{sample}/FirstSearch/RunPeptideShakerRef/PeptideShaker.log",
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/first_db_search/RunPeptideShakerRef.benchmark.txt"
    params:
        refname = "ref",
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
        "java -cp {params.peptideshaker} eu.isas.peptideshaker.cmd.PeptideShakerCLI -reference {params.refname} -fasta_file {input.ref_decoy_fasta} -identification_files {input.searchgui_zip} -spectrum_files {input.mgf} -out {output.peptide_shaker_psdb} -threads {threads} -log {log.peptide_shaker_log} > {log.stdout_log} 2> {log.stderr_log}" 


rule SimplePeptideListRef:
    input:
        peptide_shaker_psdb = RESULT_DIR / "{sample}/FirstSearch/ref.psdb"
    output:
        peptide_shaker_report = RESULT_DIR / "{sample}/FirstSearch/ref_Default_PSM_Report.txt",
    log:
        stderr_log=RESULT_DIR / "logs/{sample}/FirstSearch/SimplePeptideListRef/stderr.log",
        stdout_log=RESULT_DIR / "logs/{sample}/FirstSearch/SimplePeptideListRef/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/first_db_search/SimplePeptideListRef.benchmark.txt"
    params:
        out_dir = str(RESULT_DIR / "{sample}/FirstSearch"),
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