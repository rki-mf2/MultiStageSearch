rule copy_config:
    input:
        orig_config = "config/config.yaml"
    output:
        copied_config = RESULT_DIR / "{sample}/config.yaml",
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/utils/copy_config.benchmark.txt"
    threads: 1
    shell:
        "cp {input.orig_config} {output.copied_config}"


rule createMS2RescoreConfig:
    output:
        ms2_config = RESULT_DIR / "{sample}/MS2Rescore/ms2_config.json",
    log:
        stderr_log=RESULT_DIR / "logs/{sample}/MS2Rescore/Config/stderr.log",
        stdout_log=RESULT_DIR / "logs/{sample}/MS2Rescore/Config/stderr.log",
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/MS2Rescore/createMS2RescoreConfig.benchmark.txt"
    params:
        fragmentation_method = config["MS2Rescore"]["fragmentation_method"],
        fragment_tolerance = config["MS2Rescore"]["fragment_tolerance"],
        spectrum_pattern = config["MS2Rescore"]["spectrum_pattern"],
        threads = int(workflow.cores),
    conda:
        "../envs/base_python.yml"
    threads: 1
    script:
        "../scripts/createMS2RescoreConfig.py"