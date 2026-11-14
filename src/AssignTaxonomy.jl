module AssignTaxonomy

using FASTX, BioSequences, Kmers, Tables

include("style.jl")
include("progress.jl")
include("utils.jl")
include("classifier.jl")
include("show.jl")
export assign_taxonomy, get_reference, get_targets, names, values, ClassificationResult
end
