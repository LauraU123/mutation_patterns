LOCATIONS = ["point_mut", "one_before", "one_after", "before_after"]

rule all:
    input:
        graphs = expand("results/{location}_graph.png", location=LOCATIONS)

rule scaled_matrix:
    message:
        """Constructing scaled matrix for graph construction"""
    input:
        ref = "data/areference.gbk",
        tree = "data/rsv_a_genome.json",
    output:
        scaled_matrix = "results/{location}_matrix.csv"
    shell:
        """
        python3 scripts/constructing-matrices.py \
        --ref {input.ref} \
        --tree {input.tree} \
        --output {output.scaled_matrix}
        """

rule graph_construction:
    message:
        """Constructing graphs"""
    input:
        rules.scaled_matrix.output
    output:
        "results/{location}_graph.png"
    shell:
        """
        python3 graph.py \
        --input {input} \
        --output {output}
        """
