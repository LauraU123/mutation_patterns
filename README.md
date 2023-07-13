# RSV-A Mutation Spectrum

This workflow generates mutation spectrums for RSV-A synonymous mutations.
It includes context-based mutation spectrums.



The outputs include scaled and normalized mutation matrices and cumulative distributions of expected synonymous mutations in the RSV-A G gene duplcation.

One matrix and graph is constructed for each of the following:

* synonymous mutations in RSV-A without additional context (point_mut in the workflow)

* synonymous mutations in RSV-A based on nucleotides occurring before each mutation (one_before in the workflow)

* synonymous mutations in RSV-A based on nucleotides occurring after each mutation (one_after in the workflow)

* synonymous mutations in RSV-A based on nucleotides occurring before and after each mutation (before_after in the workflow)


### Running the workflow

The workflow can be run from the command line using Snakemake: snakemake --cores all

To change which matrices and graphs are constructed, edit the LOCATIONS list in the Snakefile.