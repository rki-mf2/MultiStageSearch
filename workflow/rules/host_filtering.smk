# from PepGM, edited

rule AddHostandCrap:
    input:
        host_fasta = HOST_FASTA,
        crap = config["contaminants"]
    output:
        concat_db = RESULT_DIR / "{sample}/Database/Host_crap.fasta"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/host_filtering/AddHostandCrap.benchmark.txt"
    shell:"cat {input} > {output}"  


rule AddDecoysHost:
    input: 
        concat_db = RESULT_DIR / "{sample}/Database/Host_crap.fasta"
    output: 
        decoy_crap_host = RESULT_DIR / "{sample}/Database/Host_crap_concatenated_target_decoy.fasta"
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/hostfiltering/HostDB/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/hostfiltering/HostDB/stdout.log" 
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/host_filtering/AddDecoysHost.benchmark.txt"
    params:
        searchgui = config["SearchGUI"],
    conda:
        "../envs/java.yml"
    threads: 1
    shell: 
        "java -cp {params.searchgui} eu.isas.searchgui.cmd.FastaCLI -in {input.concat_db} -decoy > {log.stdout_log} 2> {log.stderr_log}"


rule SearchHostSpectra:
    input: 
        mgf = MGF_FILE,
        decoy_crap_host = RESULT_DIR / "{sample}/Database/Host_crap_concatenated_target_decoy.fasta",
        par = SEARCHGUI_PAR_FILE
    output:  
        out_zip = RESULT_DIR / "{sample}/SpectraFilter/host_searchgui_out.zip"
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/hostfiltering/SearchHostSpectra/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/hostfiltering/SearchHostSpectra/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/host_filtering/SearchHostSpectra.benchmark.txt"
    params:
        searchgui = config["SearchGUI"],
        result_dir = str(RESULT_DIR / "{sample}/SpectraFilter"),
        hostname = "host",
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
        "java -cp {params.searchgui} eu.isas.searchgui.cmd.SearchCLI -spectrum_files {input.mgf} -fasta_file {input.decoy_crap_host} -output_folder {params.result_dir} -id_params {input.par} -output_default_name {params.hostname}_searchgui_out -psm_fdr {params.psm_fdr} -peptide_fdr {params.peptide_fdr} -protein_fdr {params.protein_fdr} {params.search_engine} 1 -threads {threads} > {log.stdout_log} 2> {log.stderr_log}"


rule extractSearchguiZipHost:
    input:
        searchgui_zip = RESULT_DIR / "{sample}/SpectraFilter/host_searchgui_out.zip"
    output:
        searchgui_unpacked = directory(RESULT_DIR / "{sample}/SpectraFilter/searchgui_out")
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/hostfiltering/extractSearchguiZipHost/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/hostfiltering/extractSearchguiZipHost/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/host_filtering/extractSearchguiZipHost.benchmark.txt"
    threads: 1
    shell:
        "unzip {input.searchgui_zip} -d {output.searchgui_unpacked} > {log.stdout_log}"


rule extractXTandemXMLHost:
    input:
        searchgui_unpacked = RESULT_DIR / "{sample}/SpectraFilter/searchgui_out"
    output:
        XTandem_xml =  RESULT_DIR / "{sample}/SpectraFilter/{sample}.t.xml"
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/hostfiltering/extractXTandemXMLHost/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/hostfiltering/extractXTandemXMLHost/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/host_filtering/extractXTandemXMLHost.benchmark.txt"
    params:
        XTandem_gz =  str(RESULT_DIR / "{sample}/SpectraFilter/searchgui_out/{sample}.t.xml.gz")
    threads: 1
    shell:
        "gunzip -c {params.XTandem_gz} > {output.XTandem_xml} 2> {log.stderr_log}"


rule runMS2RescoreHost:
    input:
        XTandem_xml =  RESULT_DIR / "{sample}/SpectraFilter/{sample}.t.xml",
        ms2_config = RESULT_DIR / "{sample}/MS2Rescore/ms2_config.json",
        mgf = MGF_FILE,
        decoy_crap_host = RESULT_DIR / "{sample}/Database/Host_crap_concatenated_target_decoy.fasta"
    output: 
        ms2_tsv = RESULT_DIR / "{sample}/SpectraFilter/MS2Rescore/rescored.psms.tsv"
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/hostfiltering/runMS2RescoreHost/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/hostfiltering/runMS2RescoreHost/stdout.log",
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/host_filtering/runMS2RescoreHost.benchmark.txt"
    params:
        output_path = str(RESULT_DIR / "{sample}/SpectraFilter/MS2Rescore/rescored"),
    conda:
        "../envs/MS2Rescore.yml"
    threads: workflow.cores
    shell:
        "ms2rescore -c {input.ms2_config} -s {input.mgf} -o {params.output_path} -p {input.XTandem_xml} -f {input.decoy_crap_host} > {log.stdout_log} 2> {log.stderr_log}"


rule convertMS2RescoreTSVHost:
    input:
        ms2_tsv = RESULT_DIR / "{sample}/SpectraFilter/MS2Rescore/rescored.psms.tsv",
        XTandem_xml =  RESULT_DIR / "{sample}/SpectraFilter/{sample}.t.xml"
    output:
        converted_ms2_tsv = RESULT_DIR / "{sample}/SpectraFilter/converted_ms2rescored.tsv"
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/hostfiltering/convertMS2RescoreTSV/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/hostfiltering/convertMS2RescoreTSV/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/host_filtering/convertMS2RescoreTSV.benchmark.txt"
    params:
        fdr = config["db_search"]["MS2_params"]["fdr"]
    conda:
        "../envs/pandas.yml"
    threads: 1
    script:
        "../scripts/convertMS2rescoredTsv.py"


rule RunPeptideShakerHost:
    input:
        searchgui_zip = RESULT_DIR / "{sample}/SpectraFilter/host_searchgui_out.zip",
        mgf = MGF_FILE,
        decoy_crap_host = RESULT_DIR / "{sample}/Database/Host_crap_concatenated_target_decoy.fasta"
    output: 
        peptide_shaker_psdb = RESULT_DIR / "{sample}/SpectraFilter/host.psdb"
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/hostfiltering/RunPeptideShakerHost/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/hostfiltering/RunPeptideShakerHost/stdout.log",
        peptide_shaker_log = RESULT_DIR / "logs/{sample}/hostfiltering/RunPeptideShakerHost/PeptideShaker.log",
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/host_filtering/RunPeptideShakerHost.benchmark.txt"
    params:
        hostname = "host",
        peptideshaker = config["PeptideShaker"],
    conda:
        "../envs/java.yml"
    threads: workflow.cores
    resources:
        mem_mb = config["PeptideShaker_mem_mb"]
    retries: 3 # sometimes there are java exceptions
    shell:
        "java -cp {params.peptideshaker} eu.isas.peptideshaker.cmd.PeptideShakerCLI -reference {params.hostname} -fasta_file {input.decoy_crap_host} -identification_files {input.searchgui_zip} -spectrum_files {input.mgf} -out {output.peptide_shaker_psdb} -threads {threads} -log {log.peptide_shaker_log} > {log.stdout_log} 2> {log.stderr_log}" 


rule SimplePeptideListHost:
    input:
        peptide_shaker_psdb = RESULT_DIR / "{sample}/SpectraFilter/host.psdb"
    output: 
        peptide_shaker_report = RESULT_DIR / "{sample}/SpectraFilter/host_Default_PSM_Report.txt",
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/hostfiltering/SimplePeptideListHost/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/hostfiltering/SimplePeptideListHost/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/host_filtering/SimplePeptideListHost.benchmark.txt"
    params:
        hostname = "host",
        out_dir = str(RESULT_DIR / "{sample}/SpectraFilter"),
        peptideshaker = config["PeptideShaker"],
    conda:
        "../envs/java.yml"
    threads: workflow.cores
    resources:
        mem_mb = config["PeptideShaker_mem_mb"]
    shell:
        "java -cp {params.peptideshaker} eu.isas.peptideshaker.cmd.ReportCLI -in {input.peptide_shaker_psdb} -out_reports {params.out_dir} -reports 3 > {log.stdout_log} 2> {log.stderr_log}" 
   

rule FilterSpectra:
    input: 
        mgf = MGF_FILE,
        report = SearchDB.getReports()["host_filtering_report"]
    output: 
        filtered_mgf = RESULT_DIR / "{sample}/SpectraFilter/{sample}.mgf"
    log:
        stdout_log = RESULT_DIR / "logs/{sample}/hostfiltering/FilterSpectra/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/host_filtering/FilterSpectra.benchmark.txt"
    conda: 
        "../envs/pandas.yml"
    threads: 1
    script: 
        "../scripts/host_filtering.py"