# AssignTaxonomy.jl

Taxonomic classification of DNA sequences with the
[RDP naive Bayesian classifier](https://pubmed.ncbi.nlm.nih.gov/17586664/)
(Wang et al. 2007), as in R's `dada2::assignTaxonomy`, and exact-match
species assignment, as in `dada2::assignSpecies` and `dada2::addSpecies`.
Results match R's.

```julia
using Pkg
Pkg.add("AssignTaxonomy")
```

## Usage

`seqs` is a vector of sequences, a FASTA file, or a sequence table from
[DenoiseDNA.jl](https://github.com/EvoArt/DenoiseDNA.jl). References are the
dada2-formatted training files (see
[the DADA2 reference databases](https://benjjneb.github.io/dada2/training.html)).

```julia
using AssignTaxonomy

taxa = assignTaxonomy(seqs, "silva_nr99_v138.2_toGenus_trainset.fa.gz")
taxa = addSpecies(taxa, "silva_v138.2_assignSpecies.fa.gz")
```

The result is a `TaxonomyTable`, a Tables.jl table with a `Sequence` column
and one column per rank (`missing` where the bootstrap support is below
`minBoot=50`):

```julia
using DataFrames
DataFrame(taxa)
```

Train once and reuse the model on several sets of sequences:

```julia
model = train_taxonomy("silva_nr99_v138.2_toGenus_trainset.fa.gz")
taxa1 = assignTaxonomy(seqs1, model)
taxa2 = assignTaxonomy(seqs2, model)
```

`seed` makes the bootstrap reproducible: `assignTaxonomy(seqs, ref; seed=1)`
gives R's results after `set.seed(1)`, except for sequences whose best genera
tie exactly, which R breaks at random.

`assignSpecies(seqs, ref)` gives the genus and species of each sequence by
exact matching alone.

The older `assign_taxonomy(seq_fasta, ref_fasta)`, returning a
`ClassificationResult`, is still available.


## API reference

```@autodocs
Modules = [AssignTaxonomy]
```
