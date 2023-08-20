# RSV Mutation Spectra

This workflow generates mutation spectra for RSV synonymous mutations.
It includes context-based mutation spectra for nucleotides upstream and downstream of the mutation, as well as both. This workflow can be used for either RSV-A or RSV-B. 



The outputs include scaled and normalized mutation matrices and spectra, as well as cumulative distributions of expected synonymous and nonsynonymous mutations in the RSV G gene duplcation.


### Input

The inputs for this workflow can be generated using the main [RSV Nextstrain](https://github.com/nextstrain/rsv) workflow.
They include the root sequence json, newick tree, and a tree json file annotated with mutations.



### Output

One matrix and graph is constructed for each of the following:

* synonymous mutations without additional context (point_mut in the workflow)

* synonymous mutations based on nucleotides occurring before each mutation (one_before in the workflow)

* synonymous mutations based on nucleotides occurring after each mutation (one_after in the workflow)

* synonymous mutations in RSV-A based on nucleotides occurring before and after each mutation (before_after in the workflow)


### Running the workflow

The workflow can be run from the command line using Snakemake:
```Snakemake --cores all
```

To change which matrices and graphs are constructed, edit the LOCATIONS list in the Snakefile.
