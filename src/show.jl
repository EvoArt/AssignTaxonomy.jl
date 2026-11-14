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
function _conf_hist(io, conf)
    cnt = zeros(Int, 10)
    for c in conf
        cnt[clamp(1 + floor(Int, 10c), 1, 10)] += 1
    end
    m = maximum(cnt)
    sp = join(x == 0 ? _sty(io, "·"; color=_PAL.tail) :
              _spark(io, [Float64(x)], 0, m; color=(_, _) -> _ccolor((k - 0.5) / 10)) for (k, x) in enumerate(cnt))
    return _note(io, "0 ") * sp * _note(io, " 1")
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
