# The RDP naive Bayesian classifier as implemented in R's dada2
# (assignTaxonomy / C_assign_taxonomy2), reproducing its assignments and,
# given a seed, its bootstrap confidences.

const NBOOT = 100
const MIN_REF_LEN = 20
const MIN_TAX_LEN = 50
const UNSPECIFIED = "_DADA2_UNSPECIFIED"
const DEFAULT_LEVELS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

# R computes the log-probabilities with the C library's logf; use the same one
# so the Float32 tables are bit-identical to R's on the same platform.
const _LIBM = Sys.islinux() ? "libm.so.6" : Sys.isapple() ? "libSystem.dylib" : "ucrtbase"
_logf(x::Float32) = ccall((:logf, _LIBM), Float32, (Float32,), x)

"""
    TaxonomyModel

A trained classifier: for each "genus" (unique full lineage string in the
reference) the log-probability of observing each k-mer. Made by
[`train_taxonomy`](@ref); can be passed to [`assignTaxonomy`](@ref) in place
of the reference file.
"""
struct TaxonomyModel
    k::Int
    lineages::Vector{String}      # unique taxonomy strings
    ranks::Matrix{String}         # ngenus x nlevel
    rankid::Matrix{Int32}         # ngenus x nlevel, equal ids <=> equal names at that level
    logp::Matrix{Float32}         # ngenus x 4^k (genus varies fastest)
end
nlevels(m::TaxonomyModel) = size(m.ranks, 2)
ngenera(m::TaxonomyModel) = length(m.lineages)

@inline function _nt(c::UInt8)
    c == UInt8('A') && return 0
    c == UInt8('C') && return 1
    c == UInt8('G') && return 2
    c == UInt8('T') && return 3
    return -1
end

"Sorted 0-based k-mer codes of all valid k-mer positions (duplicates kept)."
function kmer_array(seq::AbstractString, k::Int)
    b = codeunits(seq)
    out = Int32[]
    sizehint!(out, max(length(b) - k + 1, 0))
    for i in 1:length(b)-k+1
        km = 0
        for j in i:i+k-1
            nt = _nt(b[j])
            if nt < 0
                km = -1
                break
            end
            km = 4km + nt
        end
        km >= 0 && push!(out, Int32(km))
    end
    return sort!(out)
end

# R's strsplit(x, ";"): a single trailing empty field is dropped.
_rfields(t) = (f = split(t, ';'); !isempty(f) && isempty(f[end]) ? f[1:end-1] : f)

"Parse dada2-formatted reference ids into lineage strings (R's preprocessing)."
function _parse_lineages(ids::Vector{String})
    tax = String.(strip.(ids))
    if length(tax) >= 10 && all(t -> occursin(r"FU\|re[pf]s", t), tax[1:10])   # UNITE format
        @info "UNITE fungal taxonomic reference detected."
        tax = [String(split(t, '|')[5]) for t in tax]
        tax = replace.(tax, r"[pcofg]__unidentified;" => UNSPECIFIED * ";")
        tax = replace.(tax, r";s__(\w+)_" => ";s__")
        tax = replace.(tax, r";s__sp$" => ";" * UNSPECIFIED)
    end
    occursin(';', tax[1]) || error(length(split(tax[1], r"\s")) == 3 ?
        "Incorrect reference file format for assignTaxonomy (this looks like a file formatted for assignSpecies)." :
        "Incorrect reference file format for assignTaxonomy.")
    depth = [length(_rfields(t)) for t in tax]
    td = maximum(depth)
    # as R: pad by plain concatenation (so a lineage without a trailing ';' merges fields)
    tax = [depth[i] < td ? t * repeat(UNSPECIFIED * ";", td - depth[i]) : t for (i, t) in enumerate(tax)]
    return tax, td
end

"""
    train_taxonomy(ref_fasta; k=8) -> TaxonomyModel
    train_taxonomy(ref_seqs, ref_ids; k=8) -> TaxonomyModel

Train the classifier on a dada2-formatted reference (FASTA ids are
`Kingdom;Phylum;...;Genus;`). The model can be reused across calls to
[`assignTaxonomy`](@ref).
"""
function train_taxonomy(ref_fasta::AbstractString; kwargs...)
    ids, seqs = read_fasta(ref_fasta)
    return train_taxonomy(seqs, ids; kwargs...)
end

function train_taxonomy(refs::AbstractVector{<:AbstractString}, ids::AbstractVector{<:AbstractString};
                        k::Integer=8, progress::Bool=_progress_default())
    keep = length.(refs) .>= MIN_REF_LEN
    all(keep) || @warn "Some reference sequences were too short (<$(MIN_REF_LEN)nts) and were excluded."
    refs = uppercase.(String.(refs[keep])); ids = String.(ids[keep])
    tax, td = _parse_lineages(ids)
    uniq = unique(tax)
    gidx = Dict(t => i for (i, t) in enumerate(uniq))
    ref2g = [gidx[t] for t in tax]
    ng = length(uniq)
    ranks = Matrix{String}(undef, ng, td)
    for (g, t) in enumerate(uniq)
        f = _rfields(t)
        for l in 1:td
            ranks[g, l] = l <= length(f) ? f[l] : ""
        end
    end
    rankid = Matrix{Int32}(undef, ng, td)
    for l in 1:td
        lev = Dict(x => Int32(i) for (i, x) in enumerate(unique(ranks[:, l])))
        rankid[:, l] = [lev[x] for x in ranks[:, l]]
    end

    nk = 4^k
    bar = _progress(length(refs), "K-mers"; unit=" references", enabled=progress)
    # distinct k-mers of each reference, in parallel; counted serially below
    kms = Vector{Vector{Int32}}(undef, length(refs))
    Threads.@threads for r in eachindex(refs)
        kms[r] = unique!(kmer_array(refs[r], k))
        _next!(bar)
    end
    _finish!(bar)
    counts = zeros(Float32, ng, nk)
    prior = zeros(Float32, nk)
    for (r, ks) in enumerate(kms), km in ks
        counts[ref2g[r], km+1] += 1f0
        prior[km+1] += 1f0
    end
    nplus1 = ones(Float32, ng)
    for g in ref2g
        nplus1[g] += 1f0
    end
    nref = length(refs)
    @inbounds for km in 1:nk
        prior[km] = Float32((Float64(prior[km]) + 0.5) / (1.0 + nref))
    end
    Threads.@threads for km in 1:nk
        @inbounds for g in 1:ng
            counts[g, km] = _logf((counts[g, km] + prior[km]) / nplus1[g])
        end
    end
    return TaxonomyModel(k, uniq, ranks, rankid, counts)
end

# Best genus for the k-mers `ks` (summed in order, in Float32, as R does).
# Exact ties are broken uniformly at random, as R does (R's draws are not
# reproducible, so neither are tie outcomes). Returns (genus, logp, tied).
function _best_genus!(acc::Vector{Float32}, ks::AbstractVector{Int32}, m::TaxonomyModel, rng)
    fill!(acc, 0f0)
    L = m.logp
    @inbounds for km in ks
        col = Int(km) + 1
        @simd for g in eachindex(acc)
            acc[g] += L[g, col]
        end
    end
    best = 0; maxp = -floatmax(Float32); nmax = 0
    @inbounds for g in eachindex(acc)
        p = acc[g]
        if p > maxp
            maxp = p; best = g; nmax = 1
        elseif p == maxp
            nmax += 1
            rand(rng) < 1.0 / nmax && (best = g)
        end
    end
    return best, maxp, nmax > 1
end

# The classification core (C_assign_taxonomy2): best lineage of each
# sequence (0 = too short), bootstrap counts per rank, and whether an exact
# tie between genera was met anywhere for that sequence.
function _classify(model::TaxonomyModel, sq::Vector{String}; tryRC::Bool=false, nboot::Integer=NBOOT,
                   seed=nothing, multithread::Bool=true, progress::Bool=_progress_default())
    nseq = length(sq)
    nseq == 0 && error("No seqs provided to classify.")
    k = model.k; nl = nlevels(model); ng = ngenera(model)
    # one stream of uniforms laid out as R does: sequence j starts at j*max_arraylen
    max_arraylen = maximum(s -> max(length(s) - k + 1, 0), sq)
    rng = RMersenneTwister(seed === nothing ? rand(UInt32) : seed)
    unifs = runif(rng, nseq * nboot * (max_arraylen ÷ 8))
    best = zeros(Int, nseq)
    boots = zeros(Int, nseq, nl)
    ties = falses(nseq)
    bar = _progress(nseq, "Classifying"; unit=" sequences", enabled=progress)
    tiebase = seed === nothing ? rand(UInt64) : UInt64(seed)
    function classify(j, acc)
        s = sq[j]
        if length(s) >= MIN_TAX_LEN
            tie = Random.Xoshiro(hash(j, tiebase))
            ks = kmer_array(s, k)
            g, lp, t = _best_genus!(acc, ks, model, tie)
            if tryRC
                krc = kmer_array(rc(s), k)
                length(krc) == length(ks) || error("Discrepancy between forward and RC arraylen.")
                grc, lprc, trc = _best_genus!(acc, krc, model, tie)
                if lprc > lp
                    g = grc; ks = krc
                end
                t |= trc
            end
            best[j] = g
            n = length(ks); nb = n ÷ 8
            bootks = Vector{Int32}(undef, nb)
            u = (j - 1) * max_arraylen
            for _ in 1:nboot
                for i in 1:nb
                    u += 1
                    bootks[i] = ks[trunc(Int, n * unifs[u])+1]
                end
                bg, _, bt = _best_genus!(acc, bootks, model, tie)
                t |= bt
                for l in 1:nl
                    model.rankid[bg, l] == model.rankid[g, l] || break
                    boots[j, l] += 1
                end
            end
            ties[j] = t
        end
        _next!(bar)
    end
    if multithread && Threads.nthreads() > 1
        nt = Threads.nthreads()
        accs = [Vector{Float32}(undef, ng) for _ in 1:nt]
        Threads.@threads :static for t in 1:nt
            for j in t:nt:nseq
                classify(j, accs[t])
            end
        end
    else
        acc = Vector{Float32}(undef, ng)
        foreach(j -> classify(j, acc), 1:nseq)
    end
    _finish!(bar)
    return best, boots, ties
end

"""
    TaxonomyTable

Result of [`assignTaxonomy`](@ref) and [`addSpecies`](@ref): `tax[i, l]` is
the assignment of `sequences[i]` at rank `levels[l]` (`missing` when below
the bootstrap threshold or unspecified), and `boot[i, l]` the bootstrap
support out of 100 (only with `outputBootstraps=true`, otherwise it has no
columns). It is a Tables.jl table with columns `Sequence` and one per rank.
"""
struct TaxonomyTable
    sequences::Vector{String}
    levels::Vector{String}
    tax::Matrix{Union{Missing,String}}
    boot::Matrix{Int}
end
Base.size(t::TaxonomyTable) = size(t.tax)

Tables.istable(::Type{TaxonomyTable}) = true
Tables.columnaccess(::Type{TaxonomyTable}) = true
Tables.columns(t::TaxonomyTable) = t
Tables.columnnames(t::TaxonomyTable) = vcat(:Sequence, Symbol.(t.levels))
Tables.getcolumn(t::TaxonomyTable, i::Int) = i == 1 ? t.sequences : t.tax[:, i-1]
Tables.getcolumn(t::TaxonomyTable, nm::Symbol) =
    nm === :Sequence ? t.sequences : t.tax[:, findfirst(==(String(nm)), t.levels)]
Base.propertynames(t::TaxonomyTable) = Tables.columnnames(t)
Base.getproperty(t::TaxonomyTable, nm::Symbol) =
    hasfield(TaxonomyTable, nm) ? getfield(t, nm) : Tables.getcolumn(t, nm)

"""
    assignTaxonomy(seqs, ref; minBoot=50, tryRC=false, outputBootstraps=false,
                   taxLevels=["Kingdom",...,"Species"], seed=nothing, multithread=true,
                   progress) -> TaxonomyTable

Classify sequences with the RDP naive Bayesian classifier, as R's
`dada2::assignTaxonomy`. `seqs` is a vector of sequences, a FASTA file, or a
sequence table; `ref` is a dada2-formatted reference FASTA (UNITE files are
recognised) or a [`TaxonomyModel`](@ref) from [`train_taxonomy`](@ref).
Ranks with bootstrap support below `minBoot` (out of 100) are `missing`.
With `tryRC`, the reverse complement is used when it matches better.

Passing `seed` makes the bootstrap reproducible and identical to R's after
`set.seed(seed)`, except where genera tie exactly: R breaks such ties with a
non-reproducible random draw.
"""
function assignTaxonomy(seqs, ref; minBoot::Real=50, tryRC::Bool=false, outputBootstraps::Bool=false,
                        taxLevels::AbstractVector=DEFAULT_LEVELS, seed::Union{Nothing,Integer}=nothing,
                        multithread::Bool=true, progress::Bool=_progress_default())
    model = ref isa TaxonomyModel ? ref : train_taxonomy(ref; progress=progress)
    sq = _getsequences(seqs)
    minimum(length, sq) < MIN_TAX_LEN &&
        @warn "Some sequences were shorter than $MIN_TAX_LEN nts and will not receive a taxonomic classification."
    best, boots, _ = _classify(model, sq; tryRC, seed, multithread, progress)
    nl = nlevels(model)
    length(taxLevels) >= nl || throw(ArgumentError("taxLevels has fewer names than the reference has ranks ($nl)"))
    tax = Matrix{Union{Missing,String}}(missing, length(sq), nl)
    for j in eachindex(sq)
        best[j] == 0 && continue
        for l in 1:nl
            boots[j, l] >= minBoot || break
            x = model.ranks[best[j], l]
            tax[j, l] = x == UNSPECIFIED ? missing : x
        end
    end
    return TaxonomyTable(sq, String.(taxLevels[1:nl]), tax, outputBootstraps ? boots : zeros(Int, length(sq), 0))
end
