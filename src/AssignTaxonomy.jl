module AssignTaxonomy

using FASTX, BioSequences, VectorizedKmers, Tables, PrettyTables

include("utils.jl")
include("classifier.jl")
export assign_taxonomy, get_reference, get_targets, names, values, ClassificationResult
end
