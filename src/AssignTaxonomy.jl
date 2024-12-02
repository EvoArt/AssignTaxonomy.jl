module AssignTaxonomy

using FASTX, BioSequences, VectorizedKmers, Tables

include("utils.jl")
include("classifier.jl")
export assign_taxonomy,get_reference,get_targets
end
