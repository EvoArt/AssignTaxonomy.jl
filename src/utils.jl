
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


word_prior(n,N) = @fastmath (n+0.5)/(N+1)

function word_priors!(priors ::Array,N)
    for i in eachindex(priors)
        priors[i] = word_prior(priors[i],N)
    end
end

conditional_prob(m,P,M) = @fastmath (m+P)/(M+1)

#### Count kmers

function count_kmers!(counts, seq::NucSeq, ::Val{K}) where K
    @boundscheck checkbounds(counts, 1:4^K)
    K < 33 || error("Currently supports up to K=32 only")
    @inbounds for kmer in UnambiguousDNAMers{K}(seq)
        counts[kmer[1].data[1] + 1] = 1
    end
    counts;
end

count_kmers(seq::NucSeq, K) = count_kmers!(falses(4^K), seq, Val(K))

# The k-mers of each reference (in parallel), and the number of references
# containing each k-mer. The total is summed serially: adding into one shared
# array from several threads loses updates.
function count_mers(refs, k = 8; bar = nothing)
    n_refs = length(refs)
    mer_vec = Vector{BitVector}(undef, n_refs)
    Threads.@threads for i in 1:n_refs
        mer_vec[i] = count_kmers(refs[i], k)
        bar === nothing || _next!(bar)
    end
    tot_kmer_array = zeros(Float64, 4^k)
    for m in mer_vec
        tot_kmer_array .+= m
    end
    return tot_kmer_array, mer_vec
end

function count_seq_mers(seq,k = 8)
    kmer_array = count_kmers(seq,k)
    return kmer_array
end
