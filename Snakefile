LOCATIONS = ["before_after","point_mut", "one_before", "one_after"]
configfile: "config/configfile.yaml"

rule all:
    input:
        graphs = expand("results/{location}_{a_or_b}_syn_graph.png", location=LOCATIONS, a_or_b= ["a"]),
        spectra = expand("results/{location}_{a_or_b}_spectra.png", location=LOCATIONS, a_or_b= [ "a"])


rule branch_from_root:
    message:
        """Adding mutations to root sequence to reconstruct all branches"""
    input:
        root = "data/rsv_{a_or_b}_genome_root-sequence.json",
        sequences = "data/{a_or_b}_sequences.fasta",
        tree = "data/{a_or_b}_tree.nwk",
        treejson = "data/rsv_{a_or_b}_genome.json"
    output:
        reconstructed_seq = "results/{a_or_b}/reconstructed_sequences.fasta"
    shell:
        """
        python3 scripts/reconstruct_from_root.py \
        --input-root {input.root} \
        --input-tree-json {input.treejson} \
        --sequences {input.sequences} \
        --input-tree {input.tree} \
        --output {output.reconstructed_seq} 
        """

rule scaled_matrix:
    message:
        """Constructing scaled matrix for graph construction"""
    input:
        ref = "data/{a_or_b}reference.gbk",
        reconstructed = rules.branch_from_root.output,
        tree = "data/rsv_{a_or_b}_genome.json"
    params:
        type_ = lambda wildcards: f'{wildcards.location}',
        subtype = lambda wildcards: f'{wildcards.a_or_b}'
    output:
        scaled_matrix = "results/{location}_{a_or_b}_matrix.csv",
    shell:
        """
        python3 scripts/constructing-matrices.py \
        --ref {input.ref} \
        --tree {input.tree} \
        --reconstructedseq {input.reconstructed} \
        --output {output.scaled_matrix} \
        --type {params.type_} \
        --rsvsubtype {params.subtype}
        """


rule spectra:
    message:
        """Constructing mutation spectra"""
    input:
        matrix = rules.scaled_matrix.output,
    output:
        "results/{location}_{a_or_b}_spectra.png"
    params:
        type_ = lambda wildcards: f'{wildcards.location}',
        subtype = lambda wildcards: f'{wildcards.a_or_b}'
    shell:
        """
        python3 scripts/mutation_spectra_sorted.py \
        --matrix {params.type_} \
        --output {output} \
        --rsvsubtype {params.subtype}
        """


rule graph_construction:
    message:
        """Constructing graphs"""
    input:
        matrix = rules.scaled_matrix.output,
        ref = rules.scaled_matrix.input.ref,
        duplicationseq = "data/last_reconstruction_{a_or_b}.fasta",
        tree = "data/{a_or_b}_tree_duplication.nwk"
    output:
        syn = "results/{location}_{a_or_b}_syn_graph.png",
        nonsyn = "results/{location}_{a_or_b}_nonsyn_graph.png"
    params:
        type_ = lambda w: f'{w.location}',
        a_or_b = lambda w: config["reconstruct"].get(w.a_or_b)
    shell:
        """
        python3 scripts/graph.py \
        --matrix {input.matrix} \
        --syn {output.syn} \
        --nonsyn {output.nonsyn} \
        --tree {input.tree} \
        --ref {input.ref} \
        --subtype {params.a_or_b} \
        --duplicationseq {input.duplicationseq} \
        --type {params.type_}
        """
