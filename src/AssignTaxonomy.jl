module AssignTaxonomy

using FASTX, BioSequences, Tables
using CodecZlib: GzipDecompressor
import Random

include("style.jl")
include("progress.jl")
include("fasta.jl")
include("rrng.jl")
include("rdp.jl")
include("species.jl")
include("utils.jl")
include("classifier.jl")
include("show.jl")

export assignTaxonomy, assignSpecies, addSpecies, train_taxonomy, TaxonomyModel, TaxonomyTable,
       assign_taxonomy, get_reference, get_targets, ClassificationResult
end
