rule createPlots:
    input:
        Plots.get_input_Plots()
    output:
        normal_search_plots = touch(RESULT_DIR / "{sample}/Plots/create_plots.txt")
    log:
        stderr_log=RESULT_DIR / "logs/{sample}/Plots/createPlots/stderr.log",
        stdout_log=RESULT_DIR / "logs/{sample}/Plots/createPlots/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/plots/createPlots.benchmark.txt"
    params:
        sample_res_path = str(RESULT_DIR / "{sample}/Plots/"),
        sample_name = "{sample}",
        host_filtering = HOST_FILTERING,
        extra_search = EXTRA_SEARCH,
        compute_db_suitability = COMPUTE_DB_SUITABILITY
    conda:
        "../envs/plots.yml"
    threads: 1
    script:
        "../scripts/create_plots.py"


rule SimilarityHeatmap:
    input:
        orf_scores = RESULT_DIR / "{sample}/taxids/ORF_scores.csv",
        proteomes_fasta = RESULT_DIR / "{sample}/Database/filtered_proteomes.fasta",
    output:
        heatmap_plot = RESULT_DIR / "{sample}/Plots/{sample}_similarity_heatmap.png"
    log:
        stderr_log=RESULT_DIR / "logs/{sample}/Plots/SimilarityHeatmap/stderr.log",
        stdout_log=RESULT_DIR / "logs/{sample}/Plots/SimilarityHeatmap/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/plots/SimilarityHeatmap.benchmark.txt"
    params:
        num_threads = workflow.cores
    conda:
        "../envs/similarity_heatmap.yml"
    threads: workflow.cores
    script:
        "../scripts/createSimilarityHeatmap.py"


rule PeptidomeHeatmap:
    input:
        report = Plots.get_input_PeptidomeHeatmap()["report"]
    output:
        peptidome_heatmap = RESULT_DIR / "{sample}/Plots/{sample}_peptidome_heatmap.png"
    log:
        stderr_log=RESULT_DIR / "logs/{sample}/Plots/PeptidomeHeatmap/stderr.log",
        stdout_log=RESULT_DIR / "logs/{sample}/Plots/PeptidomeHeatmap/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/plots/PeptidomeHeatmap.benchmark.txt"
    conda:
        "../envs/similarity_heatmap.yml"
    threads: 1
    script:
        "../scripts/createPeptidomeHeatmap.py"


rule CreateReport:
    input:
        plots_done = RESULT_DIR / "{sample}/Plots/create_plots.txt",
        script = Plots.get_current_working_directory() +  "/workflow/scripts/create_report.Rmd",
    output:
        html_report = Plots.get_current_working_directory() / RESULT_DIR / "{sample}/ResultsReport.html"
    log:
        stderr_log=RESULT_DIR / "logs/{sample}/Plots/CreateReport/stderr.log",
        stdout_log=RESULT_DIR / "logs/{sample}/Plots/CreateReport/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/plots/CreateReport.benchmark.txt"
    params:
        sample_name = "{sample}",
        res_dir =  str(Plots.get_current_working_directory() / RESULT_DIR),
    conda:
        "../envs/report.yml"
    threads: 1
    shell:
        """
        Rscript -e "rmarkdown::render('{input.script}', output_file='{output.html_report}', params=list(sample_name='{params.sample_name}', res_dir='{params.res_dir}'))" > {log.stdout_log} 2> {log.stderr_log}
        """