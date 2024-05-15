rule CreatePepGMInput:
    input:
        PepGM.get_input_CreatePepGMInput()
    output:
        PepGM_csv = RESULT_DIR / "{sample}/PepGMInput/MultiStageSearchResults.csv",
        PepGM_input_scores = RESULT_DIR / "{sample}/PepGMInput/MultiStageSearchScores.csv",
        PepGM_strain_accessions = RESULT_DIR / "{sample}/PepGMInput/MultiStageSearchStrainAccessions.csv",
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/PepGM/CreatePepGMInput/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/PepGM/CreatePepGMInput/stdout.log",
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/pepGM/CreatePepGMInput.benchmark.txt"
    conda: 
        "../envs/pandas.yml"
    script: 
        "../scripts/create_PepGM_input.py"
    

rule CopyPhylogeneticTree:
    input:
        orig_tree = RESULT_DIR / "{sample}/Phylogeny/iqtree/phylogenetic_tree.nwk",
    output:
        copied_tree = RESULT_DIR / "{sample}/PepGMInput/phylogenetic_tree.nwk",
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/pepGM/CopyPhylogeneticTree.benchmark.txt"
    threads: 1
    shell:
        "cp {input.orig_tree} {output.copied_tree}"