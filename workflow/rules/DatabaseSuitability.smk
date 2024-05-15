rule HostFilteringNovor:
    input: 
        mgf = MGF_FILE,
        par = NOVOR_PAR_FILE
    output:  
        touch(RESULT_DIR / "{sample}/DatabaseSuitability/HostFiltering/HostFilteringNovor.done"),
        novor_dir = directory(RESULT_DIR / "{sample}/DatabaseSuitability/HostFiltering")
    log:
        stderr_log=RESULT_DIR / "logs/{sample}/DatabaseSuitability/HostFiltering/RunNovor/stderr.log",
        stdout_log=RESULT_DIR / "logs/{sample}/DatabaseSuitability/HostFiltering/RunNovor/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/DatabaseSuitability/HostFiltering/RunNovor.benchmark.txt"
    params:
        denovogui = config["DatabaseSuitability"]["DeNovoGui"],
        result_dir = str(RESULT_DIR / "{sample}/DatabaseSuitability/HostFiltering"),
    conda:
        "../envs/java.yml"
    threads: workflow.cores
    resources:
        mem_mb = config["SearchGUI_mem_mb"]
    shell: 
        "java -cp {params.denovogui} com.compomics.denovogui.cmd.DeNovoCLI -spectrum_files {input.mgf} -output_folder {params.result_dir} -id_params {input.par} -novor 1 -pepnovo 0 -directag 0 -threads {threads} > {log.stdout_log} 2> {log.stderr_log}"


rule ConcatNovorHost:
    input:
        HostFilteringNovor_finished = RESULT_DIR / "{sample}/DatabaseSuitability/HostFiltering/HostFilteringNovor.done",
        host_fasta = RESULT_DIR / "{sample}/Database/Host_crap.fasta",
    output:
        concat_fasta = RESULT_DIR / "{sample}/DatabaseSuitability/HostFiltering/Host_crap_novor.fasta"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/DatabaseSuitability/HostFiltering/ConcatNovorHost.benchmark.txt"
    params:
        novor_csv_dir = str(RESULT_DIR / "{sample}/DatabaseSuitability/HostFiltering")
    conda:
        "../envs/database_suitability.yml"
    threads: 1
    script:
        "../scripts/concat_novor.py"


rule DatabaseSuitabilityAddDecoysHost:
    input: 
        concat_db = RESULT_DIR / "{sample}/DatabaseSuitability/HostFiltering/Host_crap_novor.fasta"
    output: 
        decoy_crap_host = RESULT_DIR / "{sample}/DatabaseSuitability/HostFiltering/Host_crap_novor_concatenated_target_decoy.fasta"
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/DatabaseSuitability/HostFiltering/HostDB/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/DatabaseSuitability/HostFiltering/HostDB/stdout.log" 
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/DatabaseSuitability/HostFiltering/AddDecoysHost.benchmark.txt"
    params:
        searchgui = config["SearchGUI"],
    conda:
        "../envs/java.yml"
    threads: 1
    shell: 
        "java -cp {params.searchgui} eu.isas.searchgui.cmd.FastaCLI -in {input.concat_db} -decoy > {log.stdout_log} 2> {log.stderr_log}"


rule DatabaseSuitabilitySearchHostSpectra:
    input: 
        mgf = MGF_FILE,
        decoy_crap_host = RESULT_DIR / "{sample}/DatabaseSuitability/HostFiltering/Host_crap_novor_concatenated_target_decoy.fasta",
        par = SEARCHGUI_PAR_FILE
    output:  
        out_zip = RESULT_DIR / "{sample}/DatabaseSuitability/HostFiltering/novor_host_searchgui_out.zip"
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/DatabaseSuitability/HostFiltering/SearchHostSpectra/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/DatabaseSuitability/HostFiltering/SearchHostSpectra/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/DatabaseSuitability/HostFiltering/SearchHostSpectra.benchmark.txt"
    params:
        searchgui = config["SearchGUI"],
        result_dir = str(RESULT_DIR / "{sample}/DatabaseSuitability/HostFiltering"),
        hostname = "novor_host",
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


rule DatabaseSuitabilityextractSearchguiZipHost:
    input:
        searchgui_zip = RESULT_DIR / "{sample}/DatabaseSuitability/HostFiltering/novor_host_searchgui_out.zip"
    output:
        searchgui_unpacked = directory(RESULT_DIR / "{sample}/DatabaseSuitability/HostFiltering/novor_searchgui_out")
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/DatabaseSuitability/HostFiltering/DatabaseSuitabilityextractSearchguiZipHost/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/DatabaseSuitability/HostFiltering/DatabaseSuitabilityextractSearchguiZipHost/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/DatabaseSuitability/HostFiltering/DatabaseSuitabilityextractSearchguiZipHost.benchmark.txt"
    threads: 1
    shell:
        "unzip {input.searchgui_zip} -d {output.searchgui_unpacked} > {log.stdout_log}"


rule DatabaseSuitabilityextractXTandemXMLHost:
    input:
        searchgui_unpacked = RESULT_DIR / "{sample}/DatabaseSuitability/HostFiltering/novor_searchgui_out"
    output:
        XTandem_xml =  RESULT_DIR / "{sample}/DatabaseSuitability/HostFiltering/{sample}.t.xml"
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/DatabaseSuitability/HostFiltering/DatabaseSuitabilityextractXTandemXMLHost/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/DatabaseSuitability/HostFiltering/DatabaseSuitabilityextractXTandemXMLHost/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/DatabaseSuitability/HostFiltering/DatabaseSuitabilityextractXTandemXMLHost.benchmark.txt"
    params:
        XTandem_gz =  str(RESULT_DIR / "{sample}/DatabaseSuitability/HostFiltering/novor_searchgui_out/{sample}.t.xml.gz")
    threads: 1
    shell:
        "gunzip -c {params.XTandem_gz} > {output.XTandem_xml} 2> {log.stderr_log}"


rule DatabaseSuitabilityrunMS2RescoreHost:
    input:
        XTandem_xml =  RESULT_DIR / "{sample}/DatabaseSuitability/HostFiltering/{sample}.t.xml",
        ms2_config = RESULT_DIR / "{sample}/MS2Rescore/ms2_config.json",
        mgf = MGF_FILE,
        decoy_crap_host = RESULT_DIR / "{sample}/DatabaseSuitability/HostFiltering/Host_crap_novor_concatenated_target_decoy.fasta"
    output: 
        ms2_tsv = RESULT_DIR / "{sample}/DatabaseSuitability/HostFiltering/MS2Rescore/novor_rescored.psms.tsv"
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/DatabaseSuitability/HostFiltering/DatabaseSuitabilityrunMS2RescoreHost/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/DatabaseSuitability/HostFiltering/DatabaseSuitabilityrunMS2RescoreHost/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/DatabaseSuitability/HostFiltering/DatabaseSuitabilityrunMS2RescoreHost.benchmark.txt"
    params:
        output_path = str(RESULT_DIR / "{sample}/DatabaseSuitability/HostFiltering/MS2Rescore/novor_rescored"),
    conda:
        "../envs/MS2Rescore.yml"
    threads: workflow.cores
    shell:
        "ms2rescore -c {input.ms2_config} -s {input.mgf} -o {params.output_path} -p {input.XTandem_xml} -f {input.decoy_crap_host} > {log.stdout_log} 2> {log.stderr_log}"


rule DatabaseSuitabilityconvertMS2RescoreTSVHost:
    input:
        ms2_tsv = RESULT_DIR / "{sample}/DatabaseSuitability/HostFiltering/MS2Rescore/novor_rescored.psms.tsv",
        XTandem_xml =  RESULT_DIR / "{sample}/DatabaseSuitability/HostFiltering/{sample}.t.xml"
    output:
        converted_ms2_tsv = RESULT_DIR / "{sample}/DatabaseSuitability/HostFiltering/novor_converted_ms2rescored.tsv"
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/DatabaseSuitability/HostFiltering/DatabaseSuitabilityconvertMS2RescoreTSVHost/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/DatabaseSuitability/HostFiltering/DatabaseSuitabilityconvertMS2RescoreTSVHost/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/DatabaseSuitability/HostFiltering/DatabaseSuitabilityconvertMS2RescoreTSVHost.benchmark.txt"
    params:
        fdr = config["db_search"]["MS2_params"]["fdr"]
    conda:
        "../envs/pandas.yml"
    threads: 1
    script:
        "../scripts/convertMS2rescoredTsv.py"


rule DatabaseSuitabilityRunPeptideShakerHost:
    input:
        searchgui_zip = RESULT_DIR / "{sample}/DatabaseSuitability/HostFiltering/novor_host_searchgui_out.zip",
        mgf = MGF_FILE,
        decoy_crap_host = RESULT_DIR / "{sample}/DatabaseSuitability/HostFiltering/Host_crap_novor_concatenated_target_decoy.fasta"
    output: 
        peptide_shaker_psdb = RESULT_DIR / "{sample}/DatabaseSuitability/HostFiltering/novor_host.psdb"
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/DatabaseSuitability/HostFiltering/RunPeptideShakerHost/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/DatabaseSuitability/HostFiltering/RunPeptideShakerHost/stdout.log",
        peptide_shaker_log = RESULT_DIR / "logs/{sample}/DatabaseSuitability/HostFiltering/RunPeptideShakerHost/PeptideShaker.log",
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/DatabaseSuitability/HostFiltering/RunPeptideShakerHost.benchmark.txt"
    params:
        hostname = "novor_host",
        peptideshaker = config["PeptideShaker"],
    conda:
        "../envs/java.yml"
    threads: workflow.cores
    resources:
        mem_mb = config["PeptideShaker_mem_mb"]
    retries: 3 # sometimes there are java exceptions
    shell:
        "java -cp {params.peptideshaker} eu.isas.peptideshaker.cmd.PeptideShakerCLI -reference {params.hostname} -fasta_file {input.decoy_crap_host} -identification_files {input.searchgui_zip} -spectrum_files {input.mgf} -out {output.peptide_shaker_psdb} -threads {threads} -log {log.peptide_shaker_log} > {log.stdout_log} 2> {log.stderr_log}" 


rule DatabaseSuitabilitySimplePeptideListHost:
    input:
        peptide_shaker_psdb = RESULT_DIR / "{sample}/DatabaseSuitability/HostFiltering/novor_host.psdb"
    output: 
        peptide_shaker_report = RESULT_DIR / "{sample}/DatabaseSuitability/HostFiltering/novor_host_Default_PSM_Report.txt",
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/DatabaseSuitability/HostFiltering/SimplePeptideListHost/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/DatabaseSuitability/HostFiltering/SimplePeptideListHost/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/DatabaseSuitability/HostFiltering/SimplePeptideListHost.benchmark.txt"
    params:
        hostname = "novor_host",
        out_dir = str(RESULT_DIR / "{sample}/DatabaseSuitability/HostFiltering"),
        peptideshaker = config["PeptideShaker"],
    conda:
        "../envs/java.yml"
    threads: workflow.cores
    resources:
        mem_mb = config["PeptideShaker_mem_mb"]
    shell:
        "java -cp {params.peptideshaker} eu.isas.peptideshaker.cmd.ReportCLI -in {input.peptide_shaker_psdb} -out_reports {params.out_dir} -reports 3 > {log.stdout_log} 2> {log.stderr_log}" 
   

### ---------------------------------------------------------------------------


rule FirstSearchNovor:
    input: 
        mgf = SearchDB.get_input_MGF()["mgf"],
        par = NOVOR_PAR_FILE
    output:  
        touch(RESULT_DIR / "{sample}/DatabaseSuitability/FirstSearch/FirstSearchNovor.done"),
        novor_dir = directory(RESULT_DIR / "{sample}/DatabaseSuitability/FirstSearch")
    log:
        stderr_log=RESULT_DIR / "logs/{sample}/DatabaseSuitability/FirstSearch/RunNovor/stderr.log",
        stdout_log=RESULT_DIR / "logs/{sample}/DatabaseSuitability/FirstSearch/RunNovor/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/DatabaseSuitability/FirstSearch/RunNovor.benchmark.txt"
    params:
        denovogui = config["DatabaseSuitability"]["DeNovoGui"],
        result_dir = str(RESULT_DIR / "{sample}/DatabaseSuitability/FirstSearch"),
    conda:
        "../envs/java.yml"
    threads: workflow.cores
    resources:
        mem_mb = config["SearchGUI_mem_mb"]
    shell: 
        "java -cp {params.denovogui} com.compomics.denovogui.cmd.DeNovoCLI -spectrum_files {input.mgf} -output_folder {params.result_dir} -id_params {input.par} -novor 1 --pepnovo 0 -directag 0 -threads {threads} > {log.stdout_log} 2> {log.stderr_log}"


rule ConcatNovorRef:
    input:
        novor_done = RESULT_DIR / "{sample}/DatabaseSuitability/FirstSearch/FirstSearchNovor.done",
        ref = config["db_search"]["ref"]
    output:
        concat_fasta = RESULT_DIR / "{sample}/DatabaseSuitability/FirstSearch/novor_ref.fasta"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/DatabaseSuitability/FirstSearch/ConcatNovorRef.benchmark.txt"
    params:
        result_dir = str(RESULT_DIR / "{sample}/DatabaseSuitability/FirstSearch"),
    conda:
        "../envs/database_suitability.yml"
    threads: 1
    script:
        "../scripts/concat_novor.py"


rule DatabaseSuitabilityAddDecoysRef:
    input: 
        concat_db = RESULT_DIR / "{sample}/DatabaseSuitability/FirstSearch/novor_ref.fasta"
    output: 
        ref_decoy_fasta = RESULT_DIR / "{sample}/DatabaseSuitability/FirstSearch/novor_ref_concatenated_target_decoy.fasta"
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/DatabaseSuitability/FirstSearch/RefDB/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/DatabaseSuitability/FirstSearch/RefDB/stdout.log" 
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/DatabaseSuitability/FirstSearch/AddDecoysRef.benchmark.txt"
    params:
        searchgui = config["SearchGUI"],
    conda:
        "../envs/java.yml"
    threads: 1
    shell: 
        "java -cp {params.searchgui} eu.isas.searchgui.cmd.FastaCLI -in {input.concat_db} -decoy > {log.stdout_log} 2> {log.stderr_log}"


rule DatabaseSuitabilitySearchSpectraAgainstReference:
    input: 
        mgf = MGF_FILE,
        ref_decoy_fasta = RESULT_DIR / "{sample}/DatabaseSuitability/FirstSearch/novor_ref_concatenated_target_decoy.fasta",
        par = SEARCHGUI_PAR_FILE
    output:  
        out_zip = RESULT_DIR / "{sample}/DatabaseSuitability/FirstSearch/novor_ref_searchgui_out.zip"
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/DatabaseSuitability/FirstSearch/SearchSpectraAgainstReference/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/DatabaseSuitability/FirstSearch/SearchSpectraAgainstReference/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/DatabaseSuitability/FirstSearch/SearchSpectraAgainstReference.benchmark.txt"
    params:
        searchgui = config["SearchGUI"],
        result_dir = str(RESULT_DIR / "{sample}/DatabaseSuitability/FirstSearch"),
        refname = "novor_ref",
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
        "java -cp {params.searchgui} eu.isas.searchgui.cmd.SearchCLI -spectrum_files {input.mgf} -fasta_file {input.ref_decoy_fasta} -output_folder {params.result_dir} -id_params {input.par} -output_default_name {params.refname}_searchgui_out -psm_fdr {params.psm_fdr} -peptide_fdr {params.peptide_fdr} -protein_fdr {params.protein_fdr} {params.search_engine} 1 -threads {threads} > {log.stdout_log} 2> {log.stderr_log}"


rule DatabaseSuitabilityextractSearchguiZipFirstSearch:
    input:
        searchgui_zip = RESULT_DIR / "{sample}/DatabaseSuitability/FirstSearch/novor_ref_searchgui_out.zip"
    output:
        searchgui_unpacked = directory(RESULT_DIR / "{sample}/DatabaseSuitability/FirstSearch/novor_ref_searchgui_out")
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/DatabaseSuitability/FirstSearch/DatabaseSuitabilityextractSearchguiZipFirstSearch/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/DatabaseSuitability/FirstSearch/DatabaseSuitabilityextractSearchguiZipFirstSearch/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/DatabaseSuitability/FirstSearch/DatabaseSuitabilityextractSearchguiZipFirstSearch.benchmark.txt"
    threads: 1
    shell:
        "unzip {input.searchgui_zip} -d {output.searchgui_unpacked} > {log.stdout_log}"


rule DatabaseSuitabilityextractXTandemXMLFirstSearch:
    input:
        searchgui_unpacked = RESULT_DIR / "{sample}/DatabaseSuitability/FirstSearch/novor_ref_searchgui_out"
    output:
        XTandem_xml =  RESULT_DIR / "{sample}/DatabaseSuitability/FirstSearch/novor_{sample}.t.xml"
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/DatabaseSuitability/FirstSearch/DatabaseSuitabilityextractXTandemXMLFirstSearch/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/DatabaseSuitability/FirstSearch/DatabaseSuitabilityextractXTandemXMLFirstSearch/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/DatabaseSuitability/FirstSearch/DatabaseSuitabilityextractXTandemXMLFirstSearch.benchmark.txt"
    params:
        XTandem_gz =  str( RESULT_DIR / "{sample}/DatabaseSuitability/FirstSearch/novor_ref_searchgui_out/{sample}.t.xml.gz")
    threads: 1
    shell:
        "gunzip -c {params.XTandem_gz} > {output.XTandem_xml} 2> {log.stderr_log}"


rule DatabaseSuitabilityrunMS2RescoreFirstSearch:
    input:
        XTandem_xml =  RESULT_DIR / "{sample}/DatabaseSuitability/FirstSearch/novor_{sample}.t.xml",
        ms2_config = RESULT_DIR / "{sample}/MS2Rescore/ms2_config.json",
        mgf = MGF_FILE,
        ref_decoy_fasta = RESULT_DIR / "{sample}/DatabaseSuitability/FirstSearch/novor_ref_concatenated_target_decoy.fasta"
    output: 
        ms2_tsv = RESULT_DIR / "{sample}/DatabaseSuitability/FirstSearch/MS2Rescore/novor_rescored.psms.tsv"
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/DatabaseSuitability/FirstSearch/DatabaseSuitabilityrunMS2RescoreFirstSearch/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/DatabaseSuitability/FirstSearch/DatabaseSuitabilityrunMS2RescoreFirstSearch/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/DatabaseSuitability/FirstSearch/DatabaseSuitabilityrunMS2RescoreFirstSearch.benchmark.txt"
    params:
        output_path = str(RESULT_DIR / "{sample}/DatabaseSuitability/FirstSearch/MS2Rescore/novor_rescored"),
    conda:
        "../envs/MS2Rescore.yml"
    threads: workflow.cores
    shell:
        "ms2rescore -c {input.ms2_config} -s {input.mgf} -o {params.output_path} -p {input.XTandem_xml} -f {input.ref_decoy_fasta} > {log.stdout_log} 2> {log.stderr_log}"


rule DatabaseSuitabilityconvertMS2RescoreTSVFirstSearch:
    input:
        ms2_tsv = RESULT_DIR / "{sample}/DatabaseSuitability/FirstSearch/MS2Rescore/novor_rescored.psms.tsv",
        XTandem_xml = RESULT_DIR / "{sample}/DatabaseSuitability/FirstSearch/novor_{sample}.t.xml"
    output:
        converted_ms2_tsv = RESULT_DIR / "{sample}/DatabaseSuitability/FirstSearch/novor_converted_ms2rescored.tsv"
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/DatabaseSuitability/FirstSearch/DatabaseSuitabilityconvertMS2RescoreTSVFirstSearch/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/DatabaseSuitability/FirstSearch/DatabaseSuitabilityconvertMS2RescoreTSVFirstSearch/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/DatabaseSuitability/FirstSearch/DatabaseSuitabilityconvertMS2RescoreTSVFirstSearch.benchmark.txt"
    params:
        fdr = config["db_search"]["MS2_params"]["fdr"]
    conda:
        "../envs/pandas.yml"
    threads: 1
    script:
        "../scripts/convertMS2rescoredTsv.py"


rule DatabaseSuitabilityRunPeptideShakerRef:
    input:
        searchgui_zip = RESULT_DIR / "{sample}/DatabaseSuitability/FirstSearch/novor_ref_searchgui_out.zip",
        mgf = MGF_FILE,
        ref_decoy_fasta = RESULT_DIR / "{sample}/DatabaseSuitability/FirstSearch/novor_ref_concatenated_target_decoy.fasta"
    output: 
        peptide_shaker_psdb = RESULT_DIR / "{sample}/DatabaseSuitability/FirstSearch/novor_ref.psdb"
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/DatabaseSuitability/FirstSearch/RunPeptideShakerRef/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/DatabaseSuitability/FirstSearch/RunPeptideShakerRef/stdout.log",
        peptide_shaker_log = RESULT_DIR / "logs/{sample}/DatabaseSuitability/FirstSearch/RunPeptideShakerRef/PeptideShaker.log",
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/DatabaseSuitability/FirstSearch/RunPeptideShakerRef.benchmark.txt"
    params:
        refname = "novor_ref",
        peptideshaker = config["PeptideShaker"],
    conda:
        "../envs/java.yml"
    threads: workflow.cores
    resources:
        mem_mb = config["PeptideShaker_mem_mb"]
    retries: 3 # sometimes there are java exceptions
    shell:
        "java -cp {params.peptideshaker} eu.isas.peptideshaker.cmd.PeptideShakerCLI -reference {params.refname} -fasta_file {input.ref_decoy_fasta} -identification_files {input.searchgui_zip} -spectrum_files {input.mgf} -out {output.peptide_shaker_psdb} -threads {threads} -log {log.peptide_shaker_log} > {log.stdout_log} 2> {log.stderr_log}" 


rule DatabaseSuitabilitySimplePeptideListRef:
    input:
        peptide_shaker_psdb = RESULT_DIR / "{sample}/DatabaseSuitability/FirstSearch/novor_ref.psdb"
    output: 
        peptide_shaker_report = RESULT_DIR / "{sample}/DatabaseSuitability/FirstSearch/novor_ref_Default_PSM_Report.txt",
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/DatabaseSuitability/FirstSearch/SimplePeptideListRef/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/DatabaseSuitability/FirstSearch/SimplePeptideListRef/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/DatabaseSuitability/FirstSearch/SimplePeptideListRef.benchmark.txt"
    params:
        refname = "novor_ref",
        out_dir = str(RESULT_DIR / "{sample}/DatabaseSuitability/FirstSearch"),
        peptideshaker = config["PeptideShaker"],
    conda:
        "../envs/java.yml"
    threads: workflow.cores
    resources:
        mem_mb = config["PeptideShaker_mem_mb"]
    shell:
        "java -cp {params.peptideshaker} eu.isas.peptideshaker.cmd.ReportCLI -in {input.peptide_shaker_psdb} -out_reports {params.out_dir} -reports 3 > {log.stdout_log} 2> {log.stderr_log}" 


# ### ---------------------------------------------------------------------------


rule FinalSearchNovor:
    input: 
        mgf = SearchDB.get_input_MGF()["mgf"],
        par = NOVOR_PAR_FILE
    output:  
        touch(RESULT_DIR / "{sample}/DatabaseSuitability/FinalSearch/FinalSearchNovor.done"),
        novor_dir = directory(RESULT_DIR / "{sample}/DatabaseSuitability/FinalSearch")
    log:
        stderr_log=RESULT_DIR / "logs/{sample}/DatabaseSuitability/FinalSearch/RunNovor/stderr.log",
        stdout_log=RESULT_DIR / "logs/{sample}/DatabaseSuitability/FinalSearch/RunNovor/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/DatabaseSuitability/FinalSearch/RunNovor.benchmark.txt"
    params:
        denovogui = config["DatabaseSuitability"]["DeNovoGui"],
        result_dir = str(RESULT_DIR / "{sample}/DatabaseSuitability/FinalSearch"),
    conda:
        "../envs/java.yml"
    threads: workflow.cores
    resources:
        mem_mb = config["SearchGUI_mem_mb"]
    shell: 
        "java -cp {params.denovogui} com.compomics.denovogui.cmd.DeNovoCLI -spectrum_files {input.mgf} -output_folder {params.result_dir} -id_params {input.par} -novor 1 --pepnovo 0 -directag 0 -threads {threads} > {log.stdout_log} 2> {log.stderr_log}"


rule ConcatNovorProteomes:
    input:
        novor_done = RESULT_DIR / "{sample}/DatabaseSuitability/FinalSearch/FinalSearchNovor.done",
        proteomes = RESULT_DIR / "{sample}/Database/filtered_proteomes.fasta"
    output:
        concat_fasta = RESULT_DIR / "{sample}/DatabaseSuitability/FinalSearch/filtered_proteomes_novor.fasta"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/DatabaseSuitability/FinalSearch/ConcatNovorProteomes.benchmark.txt"
    params:
        result_dir = str(RESULT_DIR / "{sample}/DatabaseSuitability/FinalSearch"),
    conda:
        "../envs/database_suitability.yml"
    threads: 1
    script:
        "../scripts/concat_novor.py"


rule DatabaseSuitabilityAddDecoysProteome:
    input: 
        proteomes = RESULT_DIR / "{sample}/DatabaseSuitability/FinalSearch/filtered_proteomes_novor.fasta"
    output: 
        proteomes_decoy = RESULT_DIR / "{sample}/DatabaseSuitability/FinalSearch/filtered_proteomes_novor_concatenated_target_decoy.fasta"
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/DatabaseSuitability/genome2proteome/proteomeDB/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/DatabaseSuitability/genome2proteome/proteomeDB/stdout.log" 
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/DatabaseSuitability/FinalSearch/AddDecoysProteome.benchmark.txt"
    params:
        searchgui = config["SearchGUI"],
    conda:
        "../envs/java.yml"
    threads: 1
    shell: 
        "java -cp {params.searchgui} eu.isas.searchgui.cmd.FastaCLI -in {input.proteomes} -decoy > {log.stdout_log} 2> {log.stderr_log}"


rule DatabaseSuitabilityFilterDuplicateORFs:
    input:
        proteomes = RESULT_DIR / "{sample}/DatabaseSuitability/FinalSearch/filtered_proteomes_novor_concatenated_target_decoy.fasta"
    output:
        proteomes_decoy_unique = RESULT_DIR / "{sample}/DatabaseSuitability/FinalSearch/proteomes_unique_novor_concatenated_target_decoy.fasta",
        orf_list_mapping = RESULT_DIR / "{sample}/DatabaseSuitability/FinalSearch/novor_ofr_list_id_mapping.tsv",
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/DatabaseSuitability/genome2proteome/FilterDuplicates/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/DatabaseSuitability/genome2proteome/FilterDuplicates/stdout.log" 
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/DatabaseSuitability/FinalSearch/FilterDuplicateORFs.benchmark.txt"
    conda:
        "../envs/fetch_data.yml"
    threads: 1
    script:
        "../scripts/filter_and_map_duplicate_ORFs.py"


rule DatabaseSuitabilitySearchSpectraAgainstProteomes:
    input: 
        proteomes_decoy_fasta = RESULT_DIR / "{sample}/DatabaseSuitability/FinalSearch/proteomes_unique_novor_concatenated_target_decoy.fasta",
        mgf = SearchDB.get_input_MGF()["mgf"],
        par = SEARCHGUI_PAR_FILE,
    output:  
        out_zip = RESULT_DIR / "{sample}/DatabaseSuitability/FinalSearch/novor_proteomes_searchgui_out.zip"
    log:
        stderr_log=RESULT_DIR / "logs/{sample}/DatabaseSuitability/FinalSearch/SearchSpectraAgainstProteomes/stderr.log",
        stdout_log=RESULT_DIR / "logs/{sample}/DatabaseSuitability/FinalSearch/SearchSpectraAgainstProteomes/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/DatabaseSuitability/FinalSearch/SearchSpectraAgainstProteomes.benchmark.txt"
    params:
        searchgui = config["SearchGUI"],
        result_dir = str(RESULT_DIR / "{sample}/DatabaseSuitability/FinalSearch"),
        refname = "novor_proteomes",
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


rule DatabaseSuitabilityextractSearchguiZipFinalSearch:
    input:
        searchgui_zip = RESULT_DIR / "{sample}/DatabaseSuitability/FinalSearch/novor_proteomes_searchgui_out.zip"
    output:
        searchgui_unpacked = directory(RESULT_DIR / "{sample}/DatabaseSuitability/FinalSearch/novor_proteomes_searchgui_out")
    log:
        stderr_log=RESULT_DIR / "logs/{sample}/DatabaseSuitability/FinalSearch/DatabaseSuitabilityextractSearchguiZipFinalSearch/stderr.log",
        stdout_log=RESULT_DIR / "logs/{sample}/DatabaseSuitability/FinalSearch/DatabaseSuitabilityextractSearchguiZipFinalSearch/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/DatabaseSuitability/FinalSearch/DatabaseSuitabilityextractSearchguiZipFinalSearch.benchmark.txt"
    threads: 1
    shell:
        "unzip {input.searchgui_zip} -d {output.searchgui_unpacked} > {log.stdout_log}"


rule DatabaseSuitabilityextractXTandemXMLFinalSearch:
    input:
        searchgui_unpacked = RESULT_DIR / "{sample}/DatabaseSuitability/FinalSearch/novor_proteomes_searchgui_out"
    output:
        XTandem_xml = RESULT_DIR / "{sample}/DatabaseSuitability/FinalSearch/{sample}.t.xml"
    log:
        stderr_log=RESULT_DIR / "logs/{sample}/DatabaseSuitability/FinalSearch/DatabaseSuitabilityextractXTandemXMLFinalSearch/stderr.log",
        stdout_log=RESULT_DIR / "logs/{sample}/DatabaseSuitability/FinalSearch/DatabaseSuitabilityextractXTandemXMLFinalSearch/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/DatabaseSuitability/FinalSearch/DatabaseSuitabilityextractXTandemXMLFinalSearch.benchmark.txt"
    params:
        XTandem_gz =  str(RESULT_DIR / "{sample}/DatabaseSuitability/FinalSearch/novor_proteomes_searchgui_out/{sample}.t.xml.gz")
    threads: 1
    shell:
        "gunzip -c {params.XTandem_gz} > {output.XTandem_xml} 2> {log.stderr_log}"


rule DatabaseSuitabilityrunMS2RescoreFinalSearch:
    input:
        XTandem_xml =  RESULT_DIR / "{sample}/DatabaseSuitability/FinalSearch/{sample}.t.xml",
        ms2_config = RESULT_DIR / "{sample}/MS2Rescore/ms2_config.json",
        mgf = MGF_FILE,
        proteomes_decoy_fasta = RESULT_DIR / "{sample}/DatabaseSuitability/FinalSearch/proteomes_unique_novor_concatenated_target_decoy.fasta",
    output: 
        ms2_tsv = RESULT_DIR / "{sample}/DatabaseSuitability/FinalSearch/MS2Rescore/novor_rescored.psms.tsv"
    log:
        stderr_log=RESULT_DIR / "logs/{sample}/DatabaseSuitability/FinalSearch/DatabaseSuitabilityrunMS2RescoreFinalSearch/stderr.log",
        stdout_log=RESULT_DIR / "logs/{sample}/DatabaseSuitability/FinalSearch/DatabaseSuitabilityrunMS2RescoreFinalSearch/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/DatabaseSuitability/FinalSearch/DatabaseSuitabilityrunMS2RescoreFinalSearch.benchmark.txt"
    params:
        output_path = str(RESULT_DIR / "{sample}/DatabaseSuitability/FinalSearch/MS2Rescore/novor_rescored"),
    conda:
        "../envs/MS2Rescore.yml"
    threads: workflow.cores
    shell:
        "ms2rescore -c {input.ms2_config} -s {input.mgf} -o {params.output_path} -p {input.XTandem_xml} -f {input.proteomes_decoy_fasta} > {log.stdout_log} 2> {log.stderr_log}"


rule DatabaseSuitabilityconvertMS2RescoreTSVFinalSearch:
    input:
        ms2_tsv = RESULT_DIR / "{sample}/DatabaseSuitability/FinalSearch/MS2Rescore/novor_rescored.psms.tsv",
        XTandem_xml =  RESULT_DIR / "{sample}/DatabaseSuitability/FinalSearch/{sample}.t.xml",
        orf_list_mapping = RESULT_DIR / "{sample}/DatabaseSuitability/FinalSearch/novor_ofr_list_id_mapping.tsv",
    output:
        converted_ms2_tsv = RESULT_DIR / "{sample}/DatabaseSuitability/FinalSearch/novor_converted_ms2rescored.tsv"
    log:
        stderr_log=RESULT_DIR / "logs/{sample}/DatabaseSuitability/FinalSearch/DatabaseSuitabilityconvertMS2RescoreTSVFinalSearch/stderr.log",
        stdout_log=RESULT_DIR / "logs/{sample}/DatabaseSuitability/FinalSearch/DatabaseSuitabilityconvertMS2RescoreTSVFinalSearch/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/DatabaseSuitability/FinalSearch/DatabaseSuitabilityconvertMS2RescoreTSVFinalSearch.benchmark.txt"
    params:
        fdr = config["db_search"]["MS2_params"]["fdr"]
    conda:
        "../envs/pandas.yml"
    threads: 1
    script:
        "../scripts/convertMS2rescoredTsv.py"


rule DatabaseSuitabilityRunPeptideShakerProteomes:
    input:
        searchgui_zip = RESULT_DIR / "{sample}/DatabaseSuitability/FinalSearch/novor_proteomes_searchgui_out.zip",
        mgf = SearchDB.get_input_MGF()["mgf"],
        proteomes_decoy_fasta = RESULT_DIR / "{sample}/DatabaseSuitability/FinalSearch/proteomes_unique_novor_concatenated_target_decoy.fasta"
    output: 
        peptide_shaker_psdb = RESULT_DIR / "{sample}/DatabaseSuitability/FinalSearch/novor_proteomes.psdb"
    log:
        stderr_log=RESULT_DIR / "logs/{sample}/DatabaseSuitability/FinalSearch/RunPeptideShakerProteomes/stderr.log",
        stdout_log=RESULT_DIR / "logs/{sample}/DatabaseSuitability/FinalSearch/RunPeptideShakerProteomes/stdout.log",
        peptide_shaker_log=RESULT_DIR / "logs/{sample}/DatabaseSuitability/FinalSearch/RunPeptideShakerProteomes/PeptideShaker.log",
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/DatabaseSuitability/FinalSearch/RunPeptideShakerProteomes.benchmark.txt"
    params:
        refname = "novor_proteomes",
        peptideshaker = config["PeptideShaker"],
    conda:
        "../envs/java.yml"
    threads: workflow.cores
    resources:
        mem_mb = config["PeptideShaker_mem_mb"]
    retries: 3 # sometimes there are java exceptions
    shell:
        "java -cp {params.peptideshaker} eu.isas.peptideshaker.cmd.PeptideShakerCLI -reference {params.refname} -fasta_file {input.proteomes_decoy_fasta} -identification_files {input.searchgui_zip} -spectrum_files {input.mgf} -out {output.peptide_shaker_psdb} -threads {threads} -log {log.peptide_shaker_log} > {log.stdout_log} 2> {log.stderr_log}" 


rule DatabaseSuitabilitySimplePeptideListProteomes:
    input:
        peptide_shaker_psdb = RESULT_DIR / "{sample}/DatabaseSuitability/FinalSearch/novor_proteomes.psdb"
    output: 
        peptide_shaker_report = RESULT_DIR / "{sample}/DatabaseSuitability/FinalSearch/novor_proteomes_Default_PSM_Report.txt",
    log:
        stderr_log=RESULT_DIR / "logs/{sample}/DatabaseSuitability/FinalSearch/SimplePeptideListProteomes/stderr.log",
        stdout_log=RESULT_DIR / "logs/{sample}/DatabaseSuitability/FinalSearch/SimplePeptideListProteomes/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/DatabaseSuitability/FinalSearch/SimplePeptideListProteomes.benchmark.txt"
    params:
        out_dir = str(RESULT_DIR / "{sample}/DatabaseSuitability/FinalSearch"),
        peptideshaker = config["PeptideShaker"],
    conda:
        "../envs/java.yml"
    threads: workflow.cores
    resources:
        mem_mb = config["PeptideShaker_mem_mb"]
    shell:
        "java -cp {params.peptideshaker} eu.isas.peptideshaker.cmd.ReportCLI -in {input.peptide_shaker_psdb} -out_reports {params.out_dir} -reports 3 > {log.stdout_log} 2> {log.stderr_log}" 


rule DatabaseSuitabilityMapORFListsProteomes:
    input:
        peptide_shaker_report = RESULT_DIR / "{sample}/DatabaseSuitability/FinalSearch/novor_proteomes_Default_PSM_Report.txt",
        orf_list_mapping = RESULT_DIR / "{sample}/DatabaseSuitability/FinalSearch/novor_ofr_list_id_mapping.tsv",
    output: 
        converted_peptide_shaker_report = RESULT_DIR / "{sample}/DatabaseSuitability/FinalSearch/novor_converted_proteomes_Default_PSM_Report.txt",
    log:
        stderr_log=RESULT_DIR / "logs/{sample}/DatabaseSuitability/FinalSearch/DatabaseSuitabilityMapORFListsProteomes/stderr.log",
        stdout_log=RESULT_DIR / "logs/{sample}/DatabaseSuitability/FinalSearch/DatabaseSuitabilityMapORFListsProteomes/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/DatabaseSuitability/FinalSearch/DatabaseSuitabilityMapORFListsProteomes.benchmark.txt"
    conda:
        "../envs/pandas.yml"
    threads: workflow.cores
    resources:
        mem_mb = config["PeptideShaker_mem_mb"]
    script:
        "../scripts/convert_peptideshaker_report.py"


# ### ---------------------------------------------------------------------------


rule ConcatNovorProteomesExtra:
    input:
        novor_done = RESULT_DIR / "{sample}/DatabaseSuitability/FinalSearch/FinalSearchNovor.done",
        proteomes = RESULT_DIR / "{sample}/Database/extra_search_proteomes.fasta"
    output:
        concat_fasta = RESULT_DIR / "{sample}/DatabaseSuitability/ExtraSearch/proteomes_novor.fasta"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/DatabaseSuitability/ExtraSearch/ConcatNovorProteomes.benchmark.txt"
    params:
        result_dir = str(RESULT_DIR / "{sample}/DatabaseSuitability/FinalSearch"),
    conda:
        "../envs/database_suitability.yml"
    threads: 1
    script:
        "../scripts/concat_novor.py"


rule DatabaseSuitabilityAddDecoysProteomeExtra:
    input: 
        proteomes = RESULT_DIR / "{sample}/DatabaseSuitability/ExtraSearch/proteomes_novor.fasta"
    output: 
        proteomes_decoy = RESULT_DIR / "{sample}/DatabaseSuitability/ExtraSearch/proteomes_novor_concatenated_target_decoy.fasta"
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/DatabaseSuitability/extra_genome2proteome/proteomeDB/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/DatabaseSuitability/extra_genome2proteome/proteomeDB/stdout.log" 
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/DatabaseSuitability/extra_genome2proteome/AddDecoysProteome.benchmark.txt"
    params:
        searchgui = config["SearchGUI"],
    conda:
        "../envs/java.yml"
    threads: 1
    shell: 
        "java -cp {params.searchgui} eu.isas.searchgui.cmd.FastaCLI -in {input.proteomes} -decoy > {log.stdout_log} 2> {log.stderr_log}"


rule DatabaseSuitabilityFilterDuplicateORFsExtra:
    input:
        proteomes = RESULT_DIR / "{sample}/DatabaseSuitability/ExtraSearch/proteomes_novor_concatenated_target_decoy.fasta"
    output:
        proteomes_decoy_unique = RESULT_DIR / "{sample}/DatabaseSuitability/ExtraSearch/proteomes_unique_novor_concatenated_target_decoy.fasta",
        orf_list_mapping = RESULT_DIR / "{sample}/DatabaseSuitability/ExtraSearch/novor_ofr_list_id_mapping.tsv",
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/DatabaseSuitability/extra_genome2proteome/FilterDuplicates/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/DatabaseSuitability/extra_genome2proteome/FilterDuplicates/stdout.log" 
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/DatabaseSuitability/extra_genome2proteome/FilterDuplicateORFs.benchmark.txt"
    conda:
        "../envs/fetch_data.yml"
    threads: 1
    script:
        "../scripts/filter_and_map_duplicate_ORFs.py"


rule DatabaseSuitabilitySearchSpectraAgainstProteomesExtra:
    input: 
        proteomes_decoy_fasta = RESULT_DIR / "{sample}/DatabaseSuitability/ExtraSearch/proteomes_unique_novor_concatenated_target_decoy.fasta",
        mgf = SearchDB.get_input_MGF()["mgf"],
        par = SEARCHGUI_PAR_FILE,
    output:  
        out_zip = RESULT_DIR / "{sample}/DatabaseSuitability/ExtraSearch/novor_proteomes_searchgui_out.zip"
    log:
        stderr_log=RESULT_DIR / "logs/{sample}/DatabaseSuitability/ExtraSearch/SearchSpectraAgainstProteomes/stderr.log",
        stdout_log=RESULT_DIR / "logs/{sample}/DatabaseSuitability/ExtraSearch/SearchSpectraAgainstProteomes/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/DatabaseSuitability/ExtraSearch/SearchSpectraAgainstProteomes.benchmark.txt"
    params:
        searchgui = config["SearchGUI"],
        result_dir = str(RESULT_DIR / "{sample}/DatabaseSuitability/ExtraSearch"),
        refname = "novor_proteomes",
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


rule DatabaseSuitabilityextractSearchguiZipExtraSearch:
    input:
        searchgui_zip = RESULT_DIR / "{sample}/DatabaseSuitability/ExtraSearch/novor_proteomes_searchgui_out.zip"
    output:
        searchgui_unpacked = directory(RESULT_DIR / "{sample}/DatabaseSuitability/ExtraSearch/novor_proteomes_searchgui_out")
    log:
        stderr_log=RESULT_DIR / "logs/{sample}/DatabaseSuitability/ExtraSearch/DatabaseSuitabilityextractSearchguiZipExtraSearch/stderr.log",
        stdout_log=RESULT_DIR / "logs/{sample}/DatabaseSuitability/ExtraSearch/DatabaseSuitabilityextractSearchguiZipExtraSearch/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/DatabaseSuitability/ExtraSearch/DatabaseSuitabilityextractSearchguiZipExtraSearch.benchmark.txt"
    threads: 1
    shell:
        "unzip {input.searchgui_zip} -d {output.searchgui_unpacked} > {log.stdout_log}"


rule DatabaseSuitabilityextractXTandemXMLExtraSearch:
    input:
        searchgui_unpacked = RESULT_DIR / "{sample}/DatabaseSuitability/ExtraSearch/novor_proteomes_searchgui_out"
    output:
        XTandem_xml =  RESULT_DIR / "{sample}/DatabaseSuitability/ExtraSearch/{sample}.t.xml"
    log:
        stderr_log=RESULT_DIR / "logs/{sample}/DatabaseSuitability/ExtraSearch/DatabaseSuitabilityextractXTandemXMLExtraSearch/stderr.log",
        stdout_log=RESULT_DIR / "logs/{sample}/DatabaseSuitability/ExtraSearch/DatabaseSuitabilityextractXTandemXMLExtraSearch/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/DatabaseSuitability/ExtraSearch/DatabaseSuitabilityextractXTandemXMLExtraSearch.benchmark.txt"
    params:
        XTandem_gz =  str(RESULT_DIR / "{sample}/DatabaseSuitability/ExtraSearch/novor_proteomes_searchgui_out/{sample}.t.xml.gz")
    threads: 1
    shell:
        "gunzip -c {params.XTandem_gz} > {output.XTandem_xml} 2> {log.stderr_log}"


rule DatabaseSuitabilityrunMS2RescoreExtraSearch:
    input:
        XTandem_xml =  RESULT_DIR / "{sample}/DatabaseSuitability/ExtraSearch/{sample}.t.xml",
        ms2_config = RESULT_DIR / "{sample}/MS2Rescore/ms2_config.json",
        mgf = MGF_FILE,
        proteomes_decoy_fasta = RESULT_DIR / "{sample}/DatabaseSuitability/ExtraSearch/proteomes_unique_novor_concatenated_target_decoy.fasta",
    output: 
        ms2_tsv = RESULT_DIR / "{sample}/DatabaseSuitability/ExtraSearch/MS2Rescore/novor_rescored.psms.tsv"
    log:
        stderr_log=RESULT_DIR / "logs/{sample}/DatabaseSuitability/ExtraSearch/DatabaseSuitabilityrunMS2RescoreExtraSearch/stderr.log",
        stdout_log=RESULT_DIR / "logs/{sample}/DatabaseSuitability/ExtraSearch/DatabaseSuitabilityrunMS2RescoreExtraSearch/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/DatabaseSuitability/ExtraSearch/DatabaseSuitabilityrunMS2RescoreExtraSearch.benchmark.txt"
    params:
        output_path = str(RESULT_DIR / "{sample}/DatabaseSuitability/ExtraSearch/MS2Rescore/novor_rescored"),
    conda:
        "../envs/MS2Rescore.yml"
    threads: workflow.cores
    shell:
        "ms2rescore -c {input.ms2_config} -s {input.mgf} -o {params.output_path} -p {input.XTandem_xml} -f {input.proteomes_decoy_fasta} > {log.stdout_log} 2> {log.stderr_log}"


rule DatabaseSuitabilityconvertMS2RescoreTSVExtraSearch:
    input:
        ms2_tsv = RESULT_DIR / "{sample}/DatabaseSuitability/ExtraSearch/MS2Rescore/novor_rescored.psms.tsv",
        XTandem_xml =  RESULT_DIR / "{sample}/DatabaseSuitability/ExtraSearch/{sample}.t.xml",
        orf_list_mapping = RESULT_DIR / "{sample}/DatabaseSuitability/ExtraSearch/novor_ofr_list_id_mapping.tsv",
    output:
        converted_ms2_tsv = RESULT_DIR / "{sample}/DatabaseSuitability/ExtraSearch/novor_converted_ms2rescored.tsv"
    log:
        stderr_log=RESULT_DIR / "logs/{sample}/DatabaseSuitability/ExtraSearch/DatabaseSuitabilityconvertMS2RescoreTSVExtraSearch/stderr.log",
        stdout_log=RESULT_DIR / "logs/{sample}/DatabaseSuitability/ExtraSearch/DatabaseSuitabilityconvertMS2RescoreTSVExtraSearch/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/DatabaseSuitability/ExtraSearch/DatabaseSuitabilityconvertMS2RescoreTSVExtraSearch.benchmark.txt"
    params:
        fdr = config["db_search"]["MS2_params"]["fdr"]
    conda:
        "../envs/pandas.yml"
    threads: 1
    script:
        "../scripts/convertMS2rescoredTsv.py"


rule DatabaseSuitabilityRunPeptideShakerProteomesExtra:
    input:
        searchgui_zip = RESULT_DIR / "{sample}/DatabaseSuitability/ExtraSearch/novor_proteomes_searchgui_out.zip",
        mgf = SearchDB.get_input_MGF()["mgf"],
        proteomes_decoy_fasta = RESULT_DIR / "{sample}/DatabaseSuitability/ExtraSearch/proteomes_unique_novor_concatenated_target_decoy.fasta"
    output: 
        peptide_shaker_psdb = RESULT_DIR / "{sample}/DatabaseSuitability/ExtraSearch/novor_proteomes.psdb"
    log:
        stderr_log=RESULT_DIR / "logs/{sample}/DatabaseSuitability/ExtraSearch/RunPeptideShakerProteomes/stderr.log",
        stdout_log=RESULT_DIR / "logs/{sample}/DatabaseSuitability/ExtraSearch/RunPeptideShakerProteomes/stdout.log",
        peptide_shaker_log=RESULT_DIR / "logs/{sample}/DatabaseSuitability/ExtraSearch/RunPeptideShakerProteomes/PeptideShaker.log",
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/DatabaseSuitability/ExtraSearch/RunPeptideShakerProteomes.benchmark.txt"
    params:
        refname = "novor_proteomes",
        peptideshaker = config["PeptideShaker"],
    conda:
        "../envs/java.yml"
    threads: workflow.cores
    resources:
        mem_mb = config["PeptideShaker_mem_mb"]
    retries: 3 # sometimes there are java exceptions
    shell:
        "java -cp {params.peptideshaker} eu.isas.peptideshaker.cmd.PeptideShakerCLI -reference {params.refname} -fasta_file {input.proteomes_decoy_fasta} -identification_files {input.searchgui_zip} -spectrum_files {input.mgf} -out {output.peptide_shaker_psdb} -threads {threads} -log {log.peptide_shaker_log} > {log.stdout_log} 2> {log.stderr_log}" 


rule DatabaseSuitabilitySimplePeptideListProteomesExtra:
    input:
        peptide_shaker_psdb = RESULT_DIR / "{sample}/DatabaseSuitability/ExtraSearch/novor_proteomes.psdb"
    output: 
        peptide_shaker_report = RESULT_DIR / "{sample}/DatabaseSuitability/ExtraSearch/novor_proteomes_Default_PSM_Report.txt",
    log:
        stderr_log=RESULT_DIR / "logs/{sample}/DatabaseSuitability/ExtraSearch/SimplePeptideListProteomes/stderr.log",
        stdout_log=RESULT_DIR / "logs/{sample}/DatabaseSuitability/ExtraSearch/SimplePeptideListProteomes/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/DatabaseSuitability/ExtraSearch/SimplePeptideListProteomes.benchmark.txt"
    params:
        out_dir = str(RESULT_DIR / "{sample}/DatabaseSuitability/ExtraSearch"),
        peptideshaker = config["PeptideShaker"],
    conda:
        "../envs/java.yml"
    threads: workflow.cores
    resources:
        mem_mb = config["PeptideShaker_mem_mb"]
    shell:
        "java -cp {params.peptideshaker} eu.isas.peptideshaker.cmd.ReportCLI -in {input.peptide_shaker_psdb} -out_reports {params.out_dir} -reports 3 > {log.stdout_log} 2> {log.stderr_log}" 


rule DatabaseSuitabilityMapORFListsProteomesExtra:
    input:
        peptide_shaker_report = RESULT_DIR / "{sample}/DatabaseSuitability/ExtraSearch/novor_proteomes_Default_PSM_Report.txt",
        orf_list_mapping = RESULT_DIR / "{sample}/DatabaseSuitability/ExtraSearch/novor_ofr_list_id_mapping.tsv",
    output: 
        converted_peptide_shaker_report = RESULT_DIR / "{sample}/DatabaseSuitability/ExtraSearch/novor_converted_proteomes_Default_PSM_Report.txt",
    log:
        stderr_log=RESULT_DIR / "logs/{sample}/DatabaseSuitability/ExtraSearch/DatabaseSuitabilityMapORFListsProteomes/stderr.log",
        stdout_log=RESULT_DIR / "logs/{sample}/DatabaseSuitability/ExtraSearch/DatabaseSuitabilityMapORFListsProteomes/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/DatabaseSuitability/ExtraSearch/DatabaseSuitabilityMapORFListsProteomes.benchmark.txt"
    conda:
        "../envs/pandas.yml"
    threads: workflow.cores
    resources:
        mem_mb = config["PeptideShaker_mem_mb"]
    script:
        "../scripts/convert_peptideshaker_report.py"
