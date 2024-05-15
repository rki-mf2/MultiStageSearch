rule filterGenomes:
    input: 
        strain_genomes = RESULT_DIR / "{sample}/FetchData/concat_strain_genomes.fasta",
        ORF_scores = RESULT_DIR / "{sample}/taxids/ORF_scores.csv",
    output:
        used_strain_genomes = RESULT_DIR / "{sample}/Phylogeny/used_strain_genomes.fasta"
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/Phylogeny/filterGenomes/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/Phylogeny/filterGenomes/stdout.log" 
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/phylogeny/filterGenomes.benchmark.txt"
    wildcard_constraints:
        sample="[^/]+"
    conda:
        "../envs/mapping.yml"
    threads: 1
    script:
        "../scripts/filterGenomes.py"


rule MAFFT:
    input: 
        strain_genomes = RESULT_DIR / "{sample}/Phylogeny/used_strain_genomes.fasta"
    output:
        msa = RESULT_DIR / "{sample}/Phylogeny/msa.fasta"
    log:
        stderr_log = RESULT_DIR / "logs/{sample}/Phylogeny/mafft/stderr.log",
        stdout_log = RESULT_DIR / "logs/{sample}/Phylogeny/mafft/stdout.log" 
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/phylogeny/MAFFT.benchmark.txt"
    wildcard_constraints:
        sample="[^/]+"
    conda:
        "../envs/mafft.yml"
    threads: workflow.cores
    shell:
        "mafft --thread {threads} --auto {input.strain_genomes} > {output.msa} 2> {log.stderr_log}"


rule createNewick:
    input:
        msa = RESULT_DIR / "{sample}/Phylogeny/msa.fasta"
    output:
        newick_tree = RESULT_DIR / "{sample}/Phylogeny/iqtree/msa.phy.treefile",
        newick_tree_log = RESULT_DIR / "{sample}/Phylogeny/iqtree/msa.phy.log"
    log:
        stderr_log=RESULT_DIR / "logs/{sample}/Phylogeny/iqtree/stderr.log",
        stdout_log=RESULT_DIR / "logs/{sample}/Phylogeny/iqtree/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/phylogeny/createNewick.benchmark.txt"
    params:
        prefix=str(RESULT_DIR / "{sample}/Phylogeny/iqtree/msa.phy"),
    wildcard_constraints:
        sample="[^/]+"
    conda:
        "../envs/iqtree.yml"
    threads: workflow.cores
    shell:
        "iqtree -s {input.msa} -T auto --threads-max {threads} -pre {params.prefix} > {log.stdout_log} 2> {log.stderr_log}"


rule renameNodes:
    input:
        newick_tree = RESULT_DIR / "{sample}/Phylogeny/iqtree/msa.phy.treefile",
        strain_accessions = RESULT_DIR / "{sample}/FetchData/strain_accessions.csv",
    output:
        renamed_newick_tree = RESULT_DIR / "{sample}/Phylogeny/iqtree/phylogenetic_tree.nwk",
    log:
        stderr_log=RESULT_DIR / "logs/{sample}/Phylogeny/iqtree/stderr.log",
        stdout_log=RESULT_DIR / "logs/{sample}/Phylogeny/iqtree/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/phylogeny/renameNodes.benchmark.txt"
    wildcard_constraints:
        sample="[^/]+"
    conda:
        "../envs/fetch_data.yml"
    threads: 1
    script:
        "../scripts/renameNodes.py"


rule visualizeTree:
    input:
        newick_tree = RESULT_DIR / "{sample}/Phylogeny/iqtree/phylogenetic_tree.nwk",
    output:
        tree_visualization = RESULT_DIR / "{sample}/Plots/{sample}_phylogeny.png"
    log:
        stderr_log=RESULT_DIR / "logs/{sample}/Phylogeny/visualizeTree/stderr.log",
        stdout_log=RESULT_DIR / "logs/{sample}/Phylogeny/visualizeTree/stdout.log"
    benchmark:
        RESULT_DIR/ "benchmarks/{sample}/phylogeny/visualizeTree.benchmark.txt"
    wildcard_constraints:
        sample="[^/]+"
    conda:
        "../envs/fetch_data.yml"
    threads: 1
    shell:
        "export QT_QPA_PLATFORM='offscreen' && ete3 view -t {input.newick_tree} --image {output.tree_visualization} > {log.stdout_log} 2> {log.stderr_log}"