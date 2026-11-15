# Display of a ClassificationResult, in the style of DenoiseDNA's objects
# (helpers in style.jl).
#
# `show(io, x)` is one line; `show(io, MIME"text/plain"(), x)` — what the
# REPL prints — is a summary: how confident the assignments are, the most
# common taxa, and the first rows, with as many ranks as the terminal fits.

const _RANKS = (:Kingdom, :Phylum, :Class, :Order, :Family, :Genus)
_col(cr, nm::Symbol) = view(values(cr), :, lookup(cr)[nm])
_nrows(cr) = size(values(cr), 1)

# RDP's suggested cut-off is 0.8; DADA2's default (minBoot = 50) is 0.5
_ccolor(c) = c >= 0.8 ? _PAL.good : c >= 0.5 ? _PAL.amber : _PAL.bad

# (name, count), most common first
function _tally(col)
    d = Dict{String,Int}()
    for x in col
        k = string(x)
        d[k] = get(d, k, 0) + 1
    end
    return sort!(collect(d); by=p -> (-p[2], p[1]))
end

# "Pseudomonas 8 · Achromobacter 4 · … 3 more" within `w` characters
function _top_taxa(io, col, w)
    t = _tally(col)
    parts = String[]
    used = 0
    for (k, (nm, c)) in enumerate(t)
        s = nm * _note(io, " " * _fmt(c))
        rest = k < length(t) ? 12 : 0          # room for "· … N more"
        used + _vis(s) + 3 + rest > w && !isempty(parts) && break
        push!(parts, s)
        used += _vis(s) + 3
    end
    more = length(t) - length(parts)
    return join(parts, _note(io, " · ")) * (more > 0 ? _note(io, " · … $more more") : "")
end

# histogram of the confidences over 0–1 in tenths, each bar coloured by
# where its tenth falls against the cut-offs
function _conf_hist(io, conf; top="1")
    cnt = zeros(Int, 10)
    for c in conf
        cnt[clamp(1 + floor(Int, 10c), 1, 10)] += 1
    end
    m = maximum(cnt)
    sp = join(x == 0 ? _sty(io, "·"; color=_PAL.tail) :
              _spark(io, [Float64(x)], 0, m; color=(_, _) -> _ccolor((k - 0.5) / 10)) for (k, x) in enumerate(cnt))
    return _note(io, "0 ") * sp * _note(io, " " * top)
end

Base.show(io::IO, cr::ClassificationResult) =
    print(io, "ClassificationResult(", _fmt(_nrows(cr)), " sequences, ",
          _fmt(length(unique(_col(cr, :Genus)))), " genera)")

function Base.show(io::IO, ::MIME"text/plain", cr::ClassificationResult)
    n = _nrows(cr)
    genus = _col(cr, :Genus)
    _title(io, "ClassificationResult",
           "$(_fmt(n)) sequences" * (n == 0 ? "" : _note(io, " · $(_fmt(length(unique(genus)))) genera")))
    n == 0 && return
    conf = Float64.(_col(cr, :Confidence))
    sc = sort(conf)
    hi = count(>=(0.8), conf); lo = count(<(0.5), conf)
    _kv(io, "confidence", _conf_hist(io, conf) * "  median " * _sty(io, @sprintf("%.2f", sc[cld(n, 2)]); color=_ccolor(sc[cld(n, 2)])) *
                         _note(io, " · ") * _sty(io, "$(_fmt(hi)) ≥ 0.8"; color=_PAL.good) *
                         (lo > 0 ? _note(io, " · ") * _sty(io, "$(_fmt(lo)) < 0.5"; color=_PAL.bad) : ""))
    W = _termwidth(io)
    for (lab, r) in (("genera", :Genus), ("families", :Family), ("phyla", :Phylum))
        _kv(io, lab, _top_taxa(io, _col(cr, r), W - _KEY_W - 4))
    end
    # rows: ID (when there are any), then as many ranks as fit, finest first
    ids = [string(x) for x in _col(cr, :ID)]
    hasid = any(!isempty, ids)
    k = _nshow(io, n, 10)
    cells(r) = [string(x) for x in view(_col(cr, r), 1:k)]
    clip(s, w) = length(s) > w ? first(s, w - 1) * "…" : s
    # fasta ids often share a prefix and differ at the end: keep both ends
    midclip(s, w) = length(s) > w ? first(s, (w - 1) ÷ 2) * "…" * last(s, w - 1 - (w - 1) ÷ 2) : s
    idcol = hasid ? [midclip(ids[i], 24) for i in 1:k] : String[]
    confcol = [_sty(io, @sprintf("%.2f", conf[i]); color=_ccolor(conf[i])) for i in 1:k]
    used = 2 + (hasid ? max(2, maximum(length, idcol)) + 2 : 0) + 10
    shown = Symbol[]
    for r in reverse(_RANKS)
        w = min(24, max(length(string(r)), maximum(length, cells(r))))
        (r === :Genus || used + w + 2 <= W) || break
        pushfirst!(shown, r)
        used += w + 2
    end
    _section(io, "Assignments", min(W, max(used, 40)))
    rows = [vcat(hasid ? [idcol[i]] : String[], [clip(string(_col(cr, r)[i]), 24) for r in shown], confcol[i]) for i in 1:k]
    hdr = vcat(hasid ? ["ID"] : String[], string.(shown), "conf")
    _table(io, hdr, rows; align=vcat(fill(:l, length(hdr) - 1), :r))
    _more(io, n, k, "sequences")
end

# ---------------------------------------------------------------- TaxonomyTable

# "Kingdom 22 · Phylum 22 · … · Genus 18" within `w` characters
function _fit_parts(io, parts, w)
    out = String[]; used = 0
    for (k, s) in enumerate(parts)
        used + _vis(s) + 3 + (k < length(parts) ? 12 : 0) > w && !isempty(out) && break
        push!(out, s); used += _vis(s) + 3
    end
    more = length(parts) - length(out)
    return join(out, _note(io, " · ")) * (more > 0 ? _note(io, " · … $more more") : "")
end

# plural label of a rank for the "most common" lines
const _RANK_PLURAL = Dict("Kingdom" => "kingdoms", "Phylum" => "phyla", "Class" => "classes", "Order" => "orders",
                           "Family" => "families", "Genus" => "genera", "Species" => "species")
_rank_label(r) = get(_RANK_PLURAL, r, r)

Base.show(io::IO, t::TaxonomyTable) =
    print(io, "TaxonomyTable(", _fmt(length(t.sequences)), " sequences; ", join(t.levels, ", "), ")")

function Base.show(io::IO, ::MIME"text/plain", t::TaxonomyTable)
    n, nl = size(t.tax)
    _title(io, "TaxonomyTable", "$(_fmt(n)) sequences" * _note(io, " · $nl ranks"))
    (n == 0 || nl == 0) && return
    W = _termwidth(io)
    assigned = [count(!ismissing, view(t.tax, :, l)) for l in 1:nl]
    _kv(io, "assigned", _fit_parts(io, ["$(t.levels[l]) " * _sty(io, _fmt(assigned[l]); color=
                                          _ccolor(assigned[l] / n)) for l in 1:nl], W - _KEY_W - 4))
    hasboot = size(t.boot, 2) > 0
    if hasboot
        b = Float64.(t.boot[:, size(t.boot, 2)]) ./ 100
        sb = sort(b)
        _kv(io, "bootstrap", t.levels[size(t.boot, 2)] * " " * _conf_hist(io, b; top="100") * "  median " *
                             _sty(io, string(round(Int, 100sb[cld(n, 2)])); color=_ccolor(sb[cld(n, 2)])) *
                             _note(io, " · ") * _sty(io, "$(_fmt(count(>=(0.8), b))) ≥ 80"; color=_PAL.good))
    end
    # most common names at the two finest ranks that have any
    shownranks = 0
    for l in nl:-1:1
        (assigned[l] == 0 || shownranks == 2) && continue
        col = collect(skipmissing(view(t.tax, :, l)))
        _kv(io, _rank_label(t.levels[l]), _top_taxa(io, col, W - _KEY_W - 4))
        shownranks += 1
    end
    # rows: sequence, then as many ranks as fit (finest first), then the bootstrap
    k = _nshow(io, n, 10)
    cell(i, l) = (x = t.tax[i, l]; x === missing ? _sty(io, "–"; color=_PAL.tail) :
                  length(x) > 24 ? first(x, 23) * "…" : x)
    seqw = min(24, maximum(i -> length(t.sequences[i]), 1:k))
    used = 2 + seqw + 2 + (hasboot ? 6 : 0)
    shown = Int[]
    for l in nl:-1:1
        w = max(length(t.levels[l]), maximum(i -> _vis(cell(i, l)), 1:k))
        (isempty(shown) || used + w + 2 <= W) || break
        pushfirst!(shown, l)
        used += w + 2
    end
    _section(io, "Assignments", min(W, max(used, 40)))
    rows = [vcat(_seq(io, t.sequences[i], seqw), [cell(i, l) for l in shown],
                 hasboot ? [_sty(io, string(t.boot[i, size(t.boot, 2)]); color=_ccolor(t.boot[i, end] / 100))] : String[])
            for i in 1:k]
    hdr = vcat("sequence", t.levels[shown], hasboot ? ["boot"] : String[])
    _table(io, hdr, rows; align=vcat(fill(:l, length(hdr) - hasboot), hasboot ? [:r] : Symbol[]))
    _more(io, n, k, "sequences")
end

# ---------------------------------------------------------------- TaxonomyModel

Base.show(io::IO, m::TaxonomyModel) =
    print(io, "TaxonomyModel(", _fmt(ngenera(m)), " lineages, ", nlevels(m), " ranks, k=", m.k, ")")

function Base.show(io::IO, ::MIME"text/plain", m::TaxonomyModel)
    _title(io, "TaxonomyModel", "$(_fmt(ngenera(m))) lineages" * _note(io, " · $(nlevels(m)) ranks · k = $(m.k)"))
    names = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    parts = [(l <= length(names) ? names[l] : "rank $l") * " " *
             _fmt(length(unique(filter(!=(UNSPECIFIED), view(m.ranks, :, l))))) for l in 1:nlevels(m)]
    _kv(io, "names", _fit_parts(io, parts, _termwidth(io) - _KEY_W - 4) * _note(io, "  distinct per rank"))
    _kv(io, "size", _fmt(round(Int, sizeof(m.logp) / 2^20)) * " MB" * _note(io, " of k-mer log-probabilities"))
end
