
#### Read in fasta files
"""
    get_reference(ref_fasta)

Takes a path to a reference fasta. Returns a `Vector{LongDNA{4}}` of reference sequences 
and a matrix taxonomic classifications. The reference fasta must be a DADA2-formatted reference database. 
See [here](https://benjjneb.github.io/dada2/training.html) for examples.
"""
function get_reference(ref_fasta)
    taxa = []
    ref_seqs = []
    open(FASTA.Reader, ref_fasta) do reader
        for record in reader
            if length(findall(";",identifier(record))) >4
                push!(ref_seqs,LongDNA{4}(sequence(record)))
                push!(taxa,split(identifier(record),";")[1:6])
            end
        end
    end
    taxa = rotr90(hcat(taxa...))[:,end:-1:1]
    ref_seqs,taxa
end

"""
    get_targets(seq_fasta)

Takes a path to a fasta of sequences to be classified. Returns a `Vector{LongDNA{4}}` of target sequences 
and a `Vector{String}` of sequence IDs, taken from hte fasta record identifiers.
"""
function get_targets(seq_fasta)
    ids = []
    seqs = []
    open(FASTA.Reader, seq_fasta) do reader
        for record in reader
                push!(seqs,LongDNA{4}(sequence(record)))
                push!(ids,identifier(record))
        end
    end
    seqs,ids
end
