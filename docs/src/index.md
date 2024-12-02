```@meta
CurrentModule = AssignTaxonomy
```

# AssignTaxonomy

Documentation for [AssignTaxonomy](https://github.com/EvoArt/AssignTaxonomy.jl), an all Julia implementation of the [RDP Naive Bayesian Classifier algorithm](https://pubmed.ncbi.nlm.nih.gov/17586664/) for assigning taxonomic classifications based on DNA sequences. Most users will only need to use the `assign_taxonomy` function on a pair of fasta files (one with target sequences, one with a reference database). However, additional functions are provided for reading in reference and target fasta files, for those who prefer to work with Julia data structures (e.g. vectors of DNA sequences, arrays of taxonomic classifications). 

```@index
```

```@autodocs
Modules = [AssignTaxonomy]
```
