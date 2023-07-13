LOCATIONS = ["point_mut", "one_before"]

rule all:
    input:
        graphs = expand("results/{location}_graph.png", location=LOCATIONS)

rule scaled_matrix:
    message:
        """Constructing scaled matrix for graph construction"""
    input:
        ref = "data/areference.gbk",
        tree = "data/rsv_a_genome.json",
    params:
        type_ = lambda wildcards: f'{wildcards.location}'
    output:
        scaled_matrix = "results/{location}_matrix.csv"
    shell:
        """
        python3 scripts/constructing-matrices.py \
        --ref {input.ref} \
        --tree {input.tree} \
        --output {output.scaled_matrix} \
        --type {params.type_}
        """

rule graph_construction:
    message:
        """Constructing graphs"""
    input:
        matrix = rules.scaled_matrix.output,
        ref = rules.scaled_matrix.input.ref,
        duplicationseq = "data/last_reconstruction.fasta",
        tree = "data/a_tree.nwk"
    output:
        "results/{location}_graph.png"
    params:
        type_ = lambda wildcards: f'{wildcards.location}'
    shell:
        """
        python3 scripts/graph.py \
        --matrix {input.matrix} \
        --output {output} \
        --tree {input.tree} \
        --ref {input.ref} \
        --duplicationseq {input.duplicationseq} \
        --type {params.type_}
        """
