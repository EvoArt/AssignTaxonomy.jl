# The original `assign_taxonomy` API, on top of the R-faithful classifier in
# rdp.jl. Returns a `ClassificationResult`: ID, Sequence, Kingdom..Genus and
# the genus-level bootstrap confidence (0-1). For R's interface and output
# see `assignTaxonomy`.

"""
    ClassificationResult

Table returned by [`assign_taxonomy`](@ref), with columns `ID`, `Sequence`,
`Kingdom`…`Genus`, `Confidence`. Access columns as `r.Genus`, names with
`names(r)`. Tables.jl source: `DataFrame(r)`, `CSV.write(path, r)`.
"""
struct ClassificationResult{T <: AbstractVecOrMat} <: Tables.AbstractColumns
    names::Vector{Symbol}
    lookup::Dict{Symbol, Int}
    values::T
end

function classification_result(data)
    class = Symbol.(["ID", "Sequence", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Confidence"])
    ClassificationResult(class, Dict(class .=> 1:9), data)
end
classification_result(ids, seqs, assignments) = classification_result(hcat(ids, seqs, assignments))

Tables.istable(::Type{<:ClassificationResult}) = true
Base.names(m::ClassificationResult) = getfield(m, :names)
Base.values(m::ClassificationResult) = getfield(m, :values)
lookup(m::ClassificationResult) = getfield(m, :lookup)
Tables.schema(m::ClassificationResult{T}) where {T} = Tables.Schema(names(m), fill(eltype(T), size(values(m), 2)))
Tables.columnaccess(::Type{<:ClassificationResult}) = true
Tables.columns(m::ClassificationResult) = m
Tables.getcolumn(m::ClassificationResult, ::Type{T}, col::Int, nm::Symbol) where {T} = values(m)[:, col]
Tables.getcolumn(m::ClassificationResult, nm::Symbol) = values(m)[:, lookup(m)[nm]]
Tables.getcolumn(m::ClassificationResult, i::Int) = values(m)[:, i]
Tables.columnnames(m::ClassificationResult) = names(m)

const _OLD_RANKS = 6    # Kingdom..Genus

# Classify `seqs` (any sequence objects) with `model` into the 6-rank result.
function _classification(seqs, ids, model::TaxonomyModel, n_bootstrap; progress)
    sq = [uppercase(string(s)) for s in seqs]
    best, boots, _ = _classify(model, sq; nboot=n_bootstrap, progress=progress)
    nl = nlevels(model); gl = min(_OLD_RANKS, nl)
    ranks = fill("", length(sq), _OLD_RANKS)
    for j in eachindex(sq)
        best[j] == 0 && continue
        for l in 1:min(_OLD_RANKS, nl)
            x = model.ranks[best[j], l]
            ranks[j, l] = x == UNSPECIFIED ? "" : x
        end
    end
    conf = boots[:, gl] ./ n_bootstrap
    return classification_result(hcat(Any[x for x in ids], Any[s for s in seqs], ranks, conf))
end

_model(ref_fasta::AbstractString, k, lp, progress) =
    lp isa TaxonomyModel ? lp : train_taxonomy(ref_fasta; k=k, progress=progress)
function _model(refs::AbstractVector, taxa::AbstractMatrix, k, lp, progress)
    lp isa TaxonomyModel && return lp
    ids = [join(string.(taxa[i, :]), ";") * ";" for i in axes(taxa, 1)]
    return train_taxonomy([uppercase(string(r)) for r in refs], ids; k=k, progress=progress)
end

function _assign(seqs, ids, model, n_bootstrap, keep_lp, progress)
    res = _classification(seqs, ids, model, n_bootstrap; progress)
    return keep_lp ? (res, model) : res
end

"""
    assign_taxonomy(seq_fasta, ref_fasta; k=8, n_bootstrap=100, keep_lp=false, lp=false, progress)

Use the [RDP Naive Bayesian Classifier algorithm](https://pubmed.ncbi.nlm.nih.gov/17586664/) to assign taxonomic
classifications based on DNA sequence data, exactly as R's `dada2::assignTaxonomy` does. This function takes
(the paths to) two fasta files `seq_fasta` and `ref_fasta`, containing target sequences and a reference database
respectively. It returns a `Tables.jl` compatible `ClassificationResult`, containing the target sequence IDs,
target sequences, taxonomic classifications (Kingdom to Genus, without a bootstrap cut-off) and the bootstrap
confidence of the genus.

- `seq_fasta`: Path to a fasta of sequences to be classified. A vector of sequences (and optionally ids)
  can be given instead.
- `ref_fasta`: Path to a fasta reference database. A vector of reference sequences and a matrix of
  their taxonomies (as from [`get_reference`](@ref)) can be given instead.

- `k`: Length of kmers to use.
- `n_bootstrap`: Number of bootstrap iterations to perform.
- `keep_lp`: Also return the trained model (a [`TaxonomyModel`](@ref)) if `true`.
- `lp`: A trained model to use instead of training on the reference again.
- `progress`: Show progress bars for training and classification. Defaults to on
  in interactive sessions and terminals; set it for every call with
  [`AssignTaxonomy.setprogress!`](@ref).

`ref_fasta` must be a DADA2-formatted reference database.
See [here](https://benjjneb.github.io/dada2/training.html) for examples.
For R's interface (`minBoot`, `tryRC`, ...) use [`assignTaxonomy`](@ref).
"""
function assign_taxonomy(seq_fasta::AbstractString, ref_fasta; k=8, n_bootstrap=100, keep_lp=false, lp=false,
                         progress::Bool=_progress_default())
    seqs, ids = get_targets(seq_fasta)
    assign_taxonomy(seqs, ids, ref_fasta; k, n_bootstrap, keep_lp, lp, progress)
end
assign_taxonomy(seqs::AbstractVector, ref_fasta::AbstractString; kwargs...) =
    assign_taxonomy(seqs, fill("", length(seqs)), ref_fasta; kwargs...)
function assign_taxonomy(seqs::AbstractVector, ids::AbstractVector, ref_fasta::AbstractString; k=8, n_bootstrap=100,
                         keep_lp=false, lp=false, progress::Bool=_progress_default())
    _assign(seqs, ids, _model(ref_fasta, k, lp, progress), n_bootstrap, keep_lp, progress)
end
assign_taxonomy(seq::LongDNA, ref_fasta::AbstractString; kwargs...) = assign_taxonomy([seq], [""], ref_fasta; kwargs...)
assign_taxonomy(seq::LongDNA, id, ref_fasta::AbstractString; kwargs...) = assign_taxonomy([seq], [id], ref_fasta; kwargs...)

### Alternatives for working without fastas
function assign_taxonomy(seqs::AbstractVector, ids::AbstractVector, refs::AbstractVector, taxa::AbstractMatrix;
                         k=8, n_bootstrap=100, keep_lp=false, lp=false, progress::Bool=_progress_default())
    _assign(seqs, ids, _model(refs, taxa, k, lp, progress), n_bootstrap, keep_lp, progress)
end
assign_taxonomy(seqs::AbstractVector, refs::AbstractVector, taxa::AbstractMatrix; kwargs...) =
    assign_taxonomy(seqs, fill("", length(seqs)), refs, taxa; kwargs...)
assign_taxonomy(seq, id, refs::AbstractVector, taxa::AbstractMatrix; kwargs...) =
    assign_taxonomy([seq], [id], refs, taxa; kwargs...)
assign_taxonomy(seq, refs::AbstractVector, taxa::AbstractMatrix; kwargs...) =
    assign_taxonomy([seq], [""], refs, taxa; kwargs...)
