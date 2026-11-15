# Species assignment by exact matching (port of dada2's assignSpecies,
# addSpecies, mapHits and matchGenera).

# R's strsplit(x, "\\s"): split at every whitespace character; a single
# trailing empty field is dropped.
function _rsplit_ws(s::AbstractString)
    f = split(s, r"\s")
    return !isempty(f) && isempty(f[end]) ? f[1:end-1] : f
end

# 2-bit code of the first m (<= 32) bases of an ACGT string.
function _prefix_code(s::AbstractString, m::Int)
    c = UInt64(0)
    b = codeunits(s)
    @inbounds for i in 1:m
        c = (c << 2) | UInt64(_nt(b[i]))
    end
    return c
end

# For every query, the references (or, with tryRC, their reverse complements)
# that contain it exactly, like Biostrings' vcountPDict(...) > 0.
function _exact_hits(sq::Vector{String}, refs::Vector{String}, tryRC::Bool, multithread::Bool)
    m = min(minimum(length, sq), 32)
    mask = m == 32 ? typemax(UInt64) : (UInt64(1) << 2m) - 1
    index = Dict{UInt64,Vector{Int32}}()
    for (q, s) in enumerate(sq)
        push!(get!(Vector{Int32}, index, _prefix_code(s, m)), Int32(q))
    end
    qbytes = [codeunits(s) for s in sq]
    function scan!(pairs, subject, r)
        b = codeunits(subject); n = length(b)
        code = UInt64(0); valid = 0
        @inbounds for i in 1:n
            x = _nt(b[i])
            if x < 0
                valid = 0; code = UInt64(0); continue
            end
            code = ((code << 2) | UInt64(x)) & mask
            valid += 1
            valid >= m || continue
            qs = get(index, code, nothing)
            qs === nothing && continue
            start = i - m + 1
            for q in qs
                qb = qbytes[q]; L = length(qb)
                start + L - 1 <= n || continue
                ok = true
                for t in m+1:L
                    if b[start+t-1] != qb[t]
                        ok = false; break
                    end
                end
                ok && push!(pairs, (q, Int32(r)))
            end
        end
    end
    nt = multithread ? Threads.nthreads() : 1
    found = [Tuple{Int32,Int32}[] for _ in 1:nt]
    Threads.@threads :static for t in 1:nt
        for r in t:nt:length(refs)
            scan!(found[t], refs[r], r)
            tryRC && scan!(found[t], rc(refs[r]), r)
        end
    end
    hits = [Int32[] for _ in sq]
    for f in found, (q, r) in f
        push!(hits[q], r)
    end
    return [sort!(unique!(h)) for h in hits]
end

# R sorts strings with ICU's root collation (when R has ICU, as on Linux and
# macOS builds). For ASCII that is: whitespace, then punctuation in CLDR order,
# then digits, then letters ignoring case; ties are broken by case, lower first.
const _ICU_ORDER = "\t\n\v\f\r _-,;:!?.'\"()[]{}@*/\\&#%`^+<=>|~\$0123456789abcdefghijklmnopqrstuvwxyz"
const _ICU_PRIMARY = let t = fill(0xffff, 128)
    for (w, c) in enumerate(_ICU_ORDER)
        t[Int(c)+1] = w
        isletter(c) && (t[Int(uppercase(c))+1] = w)
    end
    t
end
_icu_primary(c::Char) = isascii(c) ? _ICU_PRIMARY[Int(c)+1] : UInt16(min(0xfffe, 0x0100 + UInt32(c) % 0xfe00))
_collation_key(s::AbstractString) = ([_icu_primary(c) for c in s], [isuppercase(c) for c in s])
_rsort!(v::Vector{String}) = sort!(v; by=_collation_key)

# R's mapHits: the names of the hit references (Escherichia and Shigella
# merged), "/"-joined in sorted order if there are at most `keep` of them.
function _maphits(hitrefs, names, keep)
    h = Union{Missing,String}[names[r] for r in hitrefs]
    for i in eachindex(h)
        x = h[i]
        x !== missing && (occursin("Escherichia", x) || occursin("Shigella", x)) && (h[i] = "Escherichia/Shigella")
    end
    u = unique(h)                                   # a missing name counts here...
    length(u) <= keep || return missing
    s = _rsort!(String[x for x in u if x !== missing])  # ...but sort() drops it
    return isempty(s) ? missing : join(s, "/")
end

"""
    assignSpecies(seqs, refFasta; allowMultiple=false, tryRC=false, verbose=false,
                  multithread=true) -> TaxonomyTable

Genus and species of each sequence by exact matching against a dada2
species-assignment reference (FASTA ids `ID Genus species`), as R's
`dada2::assignSpecies`. A sequence gets a genus (species) only if all the
references containing it agree on one; with `allowMultiple`, up to that many
species (all if `true`) are returned, joined with "/". Escherichia and
Shigella are reported together as "Escherichia/Shigella". Query sequences
must be ACGT only.
"""
function assignSpecies(seqs, refFasta; allowMultiple::Union{Bool,Integer}=false, tryRC::Bool=false,
                       n::Integer=2000, verbose::Bool=false, multithread::Bool=true)
    keep = allowMultiple === true ? typemax(Int) : allowMultiple === false ? 1 : Int(allowMultiple)
    sq = _getsequences(seqs)
    ids, refs = read_fasta(refFasta)
    all(isACGT, sq) || error("Non-ACGT characters present in the query sequences.")
    if length(_rsplit_ws(ids[1])) < 3
        error(count(==(';'), ids[1]) >= 3 ?
              "Incorrect reference file format for assignSpecies (this looks like a file formatted for assignTaxonomy)." :
              "Incorrect reference file format for assignSpecies.")
    end
    field(id, i) = (f = _rsplit_ws(id); i <= length(f) ? String(f[i]) : missing)
    genus = Union{Missing,String}[field(id, 2) for id in ids]
    species = Union{Missing,String}[field(id, 3) for id in ids]
    hits = _exact_hits(sq, refs, tryRC, multithread)
    tax = Matrix{Union{Missing,String}}(missing, length(sq), 2)
    for (i, h) in enumerate(hits)
        tax[i, 1] = _maphits(h, genus, 1)
        tax[i, 2] = _maphits(h, species, keep)
    end
    verbose && println(count(!ismissing, tax[:, 2]), " out of ", length(sq), " were assigned to the species level.")
    return TaxonomyTable(sq, ["Genus", "Species"], tax, zeros(Int, length(sq), 0))
end

# R's matchGenera: is the genus from exact matching consistent with the
# classifier's genus (equal, or one of its "/"- or "_"-joined parts)?
function _matchgenera(gen_tax, gen_binom; split_glyph="/")
    (gen_tax === missing || gen_binom === missing) && return false
    (isempty(gen_tax) || isempty(gen_binom)) && return false
    # as R, the genus is used as a regular expression unescaped
    return gen_tax == gen_binom ||
           occursin(Regex("^" * gen_binom * "[ _" * split_glyph * "]"), gen_tax) ||
           occursin(Regex(split_glyph * gen_binom * "\$"), gen_tax)
end

"""
    addSpecies(taxtab, refFasta; allowMultiple=false, tryRC=false, verbose=false) -> TaxonomyTable

Add a Species column to an [`assignTaxonomy`](@ref) result by exact matching
([`assignSpecies`](@ref)), keeping a species only where its genus agrees with
the table's Genus (or last) column, as R's `dada2::addSpecies`. Bootstrap
columns, if any, are kept for the original ranks.
"""
function addSpecies(taxtab::TaxonomyTable, refFasta; allowMultiple::Union{Bool,Integer}=false,
                    tryRC::Bool=false, n::Integer=2000, verbose::Bool=false, multithread::Bool=true)
    binom = assignSpecies(taxtab.sequences, refFasta; allowMultiple, tryRC, verbose, multithread)
    levels = taxtab.levels
    gcol = something(findfirst(==("Genus"), levels), length(levels))
    gen_match = [_matchgenera(taxtab.tax[i, gcol], binom.tax[i, 1]) for i in eachindex(taxtab.sequences)]
    tax = hcat(taxtab.tax, binom.tax[:, 2])
    levels = vcat(levels, "Species")
    # R blanks the first column named "Species" (the new one unless the table had one)
    tax[.!gen_match, findfirst(==("Species"), levels)] .= missing
    verbose && println("Of which ", count(!ismissing, tax[:, end]), " had genera consistent with the input table.")
    return TaxonomyTable(taxtab.sequences, levels, tax, taxtab.boot)
end
