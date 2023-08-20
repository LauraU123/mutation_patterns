# RSV Mutation Spectra Scripts


## reconstruct_from_root.py

This script reconstructs all branches given a nwk tree, an annotated json tree with mutations and tip sequences (fasta). The output is a fasta file containing all reconstructed sequences.

## contructing_matrices.py

This function constructs scaled mutation matrices for RSV. It requires several inputs, including reference genbank file and json annotated tree file.

The outputs are a series of CSV files for each context (upstream, downstream, both or point).

## mutation_spectra_sorted.py

This function takes as input the CSV matrix files, and constructs graphs to visualize their ratios. For context-dependent graphs, they are sorted based on mutation into separate graphs. 


## graph.py

This script is used to graph expected synonymous and nonsynonymous mutations occurring in the RSV duplicated regions (see also [duplication mutation analysis](https://github.com/LauraU123/duplication). 

It takes as input aligned G gene duplications, which are obtained by cutting out the duplicated region from the output of reconstruct\_from_root and aligning the sequences (as seen in the Snakefile and duplication mutation analysis workflow).

The calculation and graphing of expected mutations at each position of the duplication includes summing the relevant fractions from the mutation spectra for each context, followed by graphing. This is done separately for each context, RSV subtype and mutation type (synonymous and nonsynonymous)



