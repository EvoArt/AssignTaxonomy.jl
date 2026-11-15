# Minimal FASTA reading (optionally gzipped, multi-line records).

"""
    read_fasta(path) -> (ids, seqs)

Read a (possibly gzipped) FASTA file. `ids` are the header lines without the
'>' (not trimmed, as ShortRead's `readFasta`); `seqs` are upper-cased.
"""
function read_fasta(path::AbstractString)
    data = read(path)
    if length(data) >= 2 && data[1] == 0x1f && data[2] == 0x8b
        data = transcode(GzipDecompressor, data)
    end
    ids = String[]; seqs = String[]
    buf = IOBuffer()
    for line in eachline(IOBuffer(data))
        line = chopsuffix(line, "\r")
        if startswith(line, '>')
            isempty(ids) || push!(seqs, uppercase(String(take!(buf))))
            push!(ids, line[2:end])
        elseif !isempty(ids)
            print(buf, strip(line))
        end
    end
    isempty(ids) || push!(seqs, uppercase(String(take!(buf))))
    return ids, seqs
end

const _COMPLEMENT = let t = collect(UInt8, 0:255)
    for (a, b) in zip("ACGTMRWSYKVHDBN-acgtmrwsykvhdbn", "TGCAKYWSRMBDHVN-tgcakywsrmbdhvn")
        t[UInt8(a)+1] = UInt8(b)
    end
    t
end

"Reverse complement of a DNA string (IUPAC aware)."
function rc(s::AbstractString)
    b = codeunits(s)
    return String([_COMPLEMENT[b[i]+1] for i in reverse(eachindex(b))])
end

isACGT(s::AbstractString) = all(c -> c in (UInt8('A'), UInt8('C'), UInt8('G'), UInt8('T')), codeunits(s))

# Sequences from what users pass: a FASTA path, one sequence, a vector of
# sequences (strings or BioSequences), or an object with a `sequences` or
# `sequence` field (e.g. DenoiseDNA's SequenceTable or DadaResult).
function _getsequences(x)
    if x isa AbstractString
        return isfile(x) ? read_fasta(x)[2] : [uppercase(x)]
    elseif x isa AbstractVector
        return [uppercase(string(s)) for s in x]
    elseif hasproperty(x, :sequences)
        return _getsequences(collect(getproperty(x, :sequences)))
    elseif hasproperty(x, :sequence)
        return _getsequences(collect(getproperty(x, :sequence)))
    end
    throw(ArgumentError("Unrecognized sequence input of type $(typeof(x))"))
end
