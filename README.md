# AssignTaxonomy

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://evoart.github.io/AssignTaxonomy/index.html)
[![Build Status](https://github.com/EvoArt/AssignTaxonomy.jl/workflows/CI/badge.svg)](https://github.com/EvoArt/AssignTaxonomy.jl/actions)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/EvoArt/AssignTaxonomy.jl?svg=true)](https://ci.appveyor.com/project/EvoArt/AssignTaxonomy-jl)
[![Coverage](https://codecov.io/gh/EvoArt/AssignTaxonomy.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/EvoArt/AssignTaxonomy.jl)


AssignTaxonomy is an all Julia implementation of the [RDP Naive Bayesian Classifier algorithm](https://pubmed.ncbi.nlm.nih.gov/17586664/) for assigning taxonomic classifications based on DNA sequences. Most users will only need to use the `assign_taxonomy` function on a pair of fasta files (one with target sequences, one with a reference database). However, additional functions are provided for reading in reference and target fasta files, for those who prefer to work with Julia data structures (e.g. vectors of DNA sequences, arrays of taxonomic classifications). 


Results can be easily converted to `DataFrame` or saved to `CSV`

```
using AssignTaxonomy, CSV, DataFrames

my_results = assign_taxonomy(targets,refs)
df = DataFrame(my_results)
CSV.write("my_results.csv",my_results)
```

You can also store and reuse log_probabilities from the classifier. Basically, training the model of your reference data once and then re using it on new target data.

```
using AssignTaxonomy

my_results,my_lp = assign_taxonomy(targets,refs,keep_lp = true)
my_new_results = get_targets(some_other_target_fasta,refs,lp = my_lp)
all_my_results = AssignTaxonomy.classification_result(vcat(values(my_results),values(my_results)))
```