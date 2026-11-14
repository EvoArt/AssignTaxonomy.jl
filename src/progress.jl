# Copied from DenoiseDNA/src/progress.jl (the two packages are used together, so
# their output should look the same): keep the copies in sync.
#
# Progress bars for the long-running steps (training the classifier on the
# reference database, classifying the targets).
#
# One line per step, redrawn in place on a terminal:
#
#   Denoising   ━━━━━━━━━━━━━━━━━━━━──────────   12/20 samples   3.2 s · ETA 2.1 s
#
# and left behind as a summary when the step ends:
#
#   Denoising   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━   20/20 samples   5.3 s   1,024 ASVs
#
# When the output is not a terminal (a log file, a notebook) only that final
# line is written. Counting is thread-safe, so a bar can be advanced from
# inside `Threads.@threads`.

const _PROGRESS = Ref{Union{Nothing,Bool}}(nothing)

"""
    AssignTaxonomy.setprogress!(enabled::Bool)

Turn progress bars on or off for every function of AssignTaxonomy. By default
they are shown in interactive sessions and when `stderr` is a terminal
(but not on CI). A `progress` keyword on a call overrides this setting.
"""
setprogress!(enabled::Bool) = (_PROGRESS[] = enabled)

function _progress_default()
    p = _PROGRESS[]
    p === nothing || return p
    get(ENV, "CI", "false") == "true" && return false
    return isinteractive() || stderr isa Base.TTY
end

const _SPINNER = ('⠋', '⠙', '⠹', '⠸', '⠼', '⠴', '⠦', '⠧', '⠇', '⠏')
const _DESC_W = 13

mutable struct _Progress
    io::IO
    enabled::Bool
    tty::Bool
    desc::String
    n::Int                  # total; 0 = unknown (a spinner and a counter)
    unit::String
    count::Int
    info::String            # shown after the timing
    t0::Float64
    tlast::Float64
    spin::Int
    lock::ReentrantLock
    done::Bool
end

"""
    _progress(n, desc; unit="", enabled=_progress_default(), io=stderr)

A progress bar over `n` steps (`n = 0`: unknown, drawn as a spinner).
Advance it with [`_next!`](@ref), set its count with [`_update!`](@ref),
and close it with [`_finish!`](@ref).
"""
function _progress(n::Integer, desc::AbstractString; unit::AbstractString="",
                   enabled::Bool=_progress_default(), io::IO=stderr)
    p = _Progress(io, enabled, io isa Base.TTY, desc, n, unit, 0, "", time(), 0.0, 0, ReentrantLock(), false)
    enabled && p.tty && _render(p)
    return p
end

"Advance `p` by `k` steps (thread-safe); `info` replaces the trailing text."
function _next!(p::_Progress, k::Integer=1; info=nothing)
    p.enabled || return p
    lock(p.lock) do
        p.count += k
        info === nothing || (p.info = info)
        _maybe_render(p)
    end
    return p
end

"Set the count of `p` (e.g. clusters found so far by a spinner)."
function _update!(p::_Progress, count::Integer; info=nothing)
    p.enabled || return p
    lock(p.lock) do
        p.count = count
        info === nothing || (p.info = info)
        _maybe_render(p)
    end
    return p
end

"Set the trailing text of `p` without advancing it."
function _info!(p::_Progress, info::AbstractString)
    p.enabled || return p
    lock(p.lock) do
        p.info = info
        _maybe_render(p)
    end
    return p
end

"Close `p`, leaving its final line (with `info`) on the screen."
function _finish!(p::_Progress; info=nothing, ok::Bool=true)
    p.enabled || return p
    lock(p.lock) do
        p.done && return
        p.done = true
        info === nothing || (p.info = info)
        _render(p; final=true, ok=ok)
        println(p.io)
    end
    return p
end

"""
    _progress_println(p, msg)

Write a line of `verbose` output above a running bar, so the two do not
end up on the same line.
"""
function _progress_println(p::Union{Nothing,_Progress}, msg)
    if p === nothing || !p.enabled || !p.tty || p.done
        println(msg)
        return
    end
    lock(p.lock) do
        print(p.io, "\r\e[2K")
        flush(p.io)
        println(msg)
        flush(stdout)
        _render(p)
    end
end

function _maybe_render(p::_Progress)
    (p.tty && !p.done) || return
    t = time()
    (t - p.tlast >= 0.1 || (p.n > 0 && p.count >= p.n)) || return
    _render(p)
end

_secs(t) = t < 1 ? "$(round(Int, 1000t)) ms" : t < 60 ? "$(round(t; digits=1)) s" :
           t < 3600 ? "$(floor(Int, t / 60)) min $(round(Int, mod(t, 60))) s" :
           "$(floor(Int, t / 3600)) h $(round(Int, mod(t, 3600) / 60)) min"

function _render(p::_Progress; final::Bool=false, ok::Bool=true)
    io = p.io
    p.tlast = time()
    el = p.tlast - p.t0
    unit = isempty(p.unit) ? "" : " " * strip(p.unit)
    desc = _sty(io, rpad(p.desc, _DESC_W); color=final ? (ok ? _PAL.good : _PAL.warn) : _PAL.head, bold=true)
    if p.n > 0
        frac = clamp(p.count / p.n, 0.0, 1.0)
        cnt = final && p.count < p.n ? "$(_fmt(p.count))$unit" : "$(_fmt(p.count))/$(_fmt(p.n))$unit"
        tim = final || p.count == 0 ? _secs(el) :
              p.count >= p.n ? _secs(el) : _secs(el) * " · ETA " * _secs(el / p.count * (p.n - p.count))
        right = "  " * lpad(cnt, 2 * ndigits(p.n) + 2 + length(unit)) * "  " * _sty(io, rpad(tim, 7); color=_PAL.key)
        info = isempty(p.info) ? "" : "  " * p.info
        # a fixed width (not what the info leaves), so finished bars line up
        bw = clamp(displaysize(io)[2] - _DESC_W - 60, 10, 24)
        full = final ? bw : floor(Int, frac * bw)
        bar = _sty(io, repeat('━', full); color=final ? (ok ? _PAL.good : _PAL.warn) : _PAL.head) *
              _sty(io, repeat('─', bw - full); color=_PAL.rule)
        line = desc * "  " * bar * right * info
    else
        p.spin = mod1(p.spin + 1, length(_SPINNER))
        mark = final ? _sty(io, ok ? "✓" : "✗"; color=ok ? _PAL.good : _PAL.warn) :
               _sty(io, string(_SPINNER[p.spin]); color=_PAL.head)
        cnt = p.count > 0 || !isempty(unit) ? "  " * _fmt(p.count) * unit : ""
        line = desc * "  " * mark * cnt * "  " * _sty(io, _secs(el); color=_PAL.key) *
               (isempty(p.info) ? "" : "  " * p.info)
    end
    if p.tty
        print(io, "\r\e[2K", rstrip(line))
        flush(io)
    elseif final
        print(io, rstrip(line))
    end
    return nothing
end
