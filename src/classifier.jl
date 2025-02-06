
#https://doi.org/10.1128/AEM.00062-07
#The Classifier uses a feature space consisting of all possible 8-base
#subsequences (words).
# The position of the word in a sequence is ignored.
# Only the words occurring in the query contribute to the score.

# Most compute time is spent on assignment (largely within the bootstrapping process). This 
# is where most effort should be spent, outside of ensuring correctness and usability.
# Memory usage is also worth looking at, to ensure that users can run the code in parallel
# across many threads without issues.

"""
The results returned by `assign_taxonomy`. Individual columns can be accessed by e.g. `my_result.Genus`, 
and a list of column names can be accessed by `names(my_result)`. The result can 
be converted to a `DataFrame` by `using DataFrames; DataFrame(my_result)` or written to CSV with headers by 
`using CSV; CSV.write("my_result.csv",my_result)`.
"""
struct ClassificationResult{T <: AbstractVecOrMat} <: Tables.AbstractColumns
    names::Vector{Symbol}
    lookup::Dict{Symbol, Int}
    values::T
end

function classification_result(data) 
    class = Symbol.(["ID","Sequence","Kingdom", "Phylum", "Class", "Order", "Family", "Genus","Confidence"])
    ClassificationResult(class,
    Dict(class .=> 1:9),
    data
    )
end

function classification_result(ids,seqs,assignments) 
    classification_result(hcat(ids,seqs,assignments))
   
end

function Base.show(io::IO, cr::ClassificationResult)
    pretty_table(
    values(cr)[:,[1,8,9]];
    formatters    = ft_printf("%5.2f", 2:4),
    header        = names(cr)[[1,8,9]],
    header_crayon = crayon"yellow bold")
end

Tables.istable(::Type{<:ClassificationResult}) = true
names(m::ClassificationResult) = getfield(m, :names)
values(m::ClassificationResult) = getfield(m, :values)
lookup(m::ClassificationResult) = getfield(m, :lookup)
Tables.schema(m::ClassificationResult{T}) where {T} = Tables.Schema(names(m), fill(eltype(T), size(values(m), 2)))
Tables.columnaccess(::Type{<:ClassificationResult}) = true
Tables.columns(m::ClassificationResult) = m
Tables.getcolumn(m::ClassificationResult, ::Type{T}, col::Int, nm::Symbol) where {T} = values(m)[:, col]
Tables.getcolumn(m::ClassificationResult, nm::Symbol) = values(m)[:, lookup(m)[nm]]
Tables.getcolumn(m::ClassificationResult, i::Int) = values(m)[:, i]
Tables.columnnames(m::ClassificationResult) = names(m)




#### Underlying algorithm
function naieve_bayes(seqs::Vector,refs::Vector,k, n_bootstrap,lp,genera)
    N = length(refs)
    n = length(seqs)
    gens = sort(unique(genera))
    g = length(gens)
    counts = zeros(g)
    gen2ind = Dict(gens .=> 1:g)
    seq2gen = Dict(1:N .=> [gen2ind[gen] for gen in genera])
    assignments = Vector{Int64}(undef,n)
    confs = Vector{Float64}(undef,n)
    if lp == false
        # This gives correct results on test cases
        # but could likely be sped up.
        log_probs = [zeros(Float32,4^k) for _ in 1:g]
        priors, a =count_mers(refs)
        word_priors!(priors,N) 
        Threads.@threads for i in 1:N
            @fastmath log_probs[seq2gen[i]] .+= conditional_prob.(a[i],priors,1)
            counts[seq2gen[i]] += 1
        end
        for i in 1:g
            @fastmath log_probs[i] .= log.(log_probs[i] ./ counts[i])
        end
    else
        log_probs = lp
    end
    Threads.@threads for i in 1:n
        kmer_array = count_seq_mers(seqs[i])
        assignment = assign(eachindex(kmer_array)[kmer_array],log_probs)
        assignments[i] =assignment
        sample_size = sum(kmer_array) ÷ k
        confs[i] = bootstrap(vec(kmer_array),log_probs,gens[assignment],sample_size,n_bootstrap,gens)
    end
    return assignments, confs, log_probs, gens
end

function naieve_bayes(seqs::Vector,refs::Vector,taxa ::Array,k, n_bootstrap,lp=false)
    # This and sbove need a bit of a tidy up, following switch to mean over genera 
    a,c,l,g = naieve_bayes(seqs,refs,k, n_bootstrap,lp,taxa[:,end])
    gens = g[a]
    inds = [findfirst(x -> x== gen,taxa[:,end]) for gen in gens]
    t = taxa[inds,:]
    return hcat(t,c),l
end

function assign(seq_mask,log_probs)
    # This dominates compute time. optimization of other parts of the code base which don't
    # speed up the assignment function are largely redundant at this point.
    cp_max = Float32(-Inf)
    ind = 1
    for i in eachindex(log_probs)
        cp = Float32(0)
        @simd for j in seq_mask
            @fastmath cp += log_probs[i][j]
        end
        if cp > cp_max
            cp_max = cp
            ind = i
        end
    end
    return ind
end


function bootstrap(kmer_vec,log_probs,assignment, sample_size,n_bootstrap,genera) 

    hits = 0
    seq_inds = eachindex(kmer_vec)[kmer_vec]
    inds =  seq_inds[1:sample_size] 
    for i in 1:n_bootstrap
        for j in 1:sample_size
            inds[j] = rand(seq_inds)
        end
        if genera[assign(inds,log_probs)]  == assignment
            hits +=1
        end
    end
    return hits/n_bootstrap
end

#### Top level function

### Most typical function call, taking reference and target fasta input

"""
    assign_taxonomy(seq_fasta,ref_fasta; k = 8, n_bootstrap = 100,keep_lp = false,lp=false)


Use the [RDP Naive Bayesian Classifier algorithm](https://pubmed.ncbi.nlm.nih.gov/17586664/) to assign taxonomic 
classifications based on DNA sequence data. This function takes (the paths to) two fasta files
`seq_fasta` and `ref_fasta`, containing target sequences and a reference database respectively.
It  returns `Tables.jl` compatible `ClassificationResult`, containing the target sequence IDs, 
target sequences, taxonomic classifications and bootstrapped confidence levels.

- `seq_fasta`: Path to a fasta of sequnces to be classified.
- `ref_fasta`: Path to a fasta reference database.

- `k`: Length of kmers to use.
- `n_bootstrap`: Number of bootstrap iterations to perform.
- `keep_lp`: Return array of log probabilities alongside classification result if `true` 
- `lp`: Array of log probabilities for the classifier to use. 

`ref_fasta` must be a DADA2-formatted reference database. 
See [here](https://benjjneb.github.io/dada2/training.html) for examples.
"""
function assign_taxonomy(seq_fasta,ref_fasta; k = 8, n_bootstrap = 100,keep_lp = false,lp=false)

    seqs,ids = get_targets(seq_fasta)
    refs,taxa = get_reference(ref_fasta)
    assign_taxonomy(seqs,ids,refs,taxa,k=k,n_bootstrap=n_bootstrap,keep_lp = keep_lp,lp = lp)
end

### Alternatives with reference fasta and LongDNA/vector targets
function assign_taxonomy(seqs::Vector,ref_fasta; k = 8, n_bootstrap = 100,keep_lp = false,lp=false)
    refs,taxa = get_reference(ref_fasta)
    assign_taxonomy(seqs,refs,taxa,k=k,n_bootstrap=n_bootstrap,keep_lp = keep_lp,lp = lp)
end
function assign_taxonomy(seqs::Vector, ids::Vector, ref_fasta; k = 8, n_bootstrap = 100,keep_lp = false,lp=false)
    refs,taxa = get_reference(ref_fasta)
    assign_taxonomy(seqs,ids,refs,taxa,k=k,n_bootstrap=n_bootstrap,keep_lp = keep_lp,lp = lp)
end
function assign_taxonomy(seq::LongDNA,ref_fasta; k = 8, n_bootstrap = 100,keep_lp = false,lp=false)
    refs,taxa = get_reference(ref_fasta)
    assign_taxonomy(seq,refs,taxa,k=k,n_bootstrap=n_bootstrap,keep_lp = keep_lp,lp = lp)
end
function assign_taxonomy(seq::LongDNA,id,ref_fasta; k = 8, n_bootstrap = 100,keep_lp = false,lp=false)
    refs,taxa = get_reference(ref_fasta)
    assign_taxonomy(seq,id,refs,taxa,k=k,n_bootstrap=n_bootstrap,keep_lp = keep_lp,lp = lp)
end

### Alternatives for working without fastas
function assign_taxonomy(seqs::Vector,ids,refs,taxa; k = 8, n_bootstrap = 100,keep_lp = false,lp=false)
    assignments,log_probs = naieve_bayes(seqs,refs,taxa,k,n_bootstrap,lp)
    res = classification_result(ids,seqs,assignments)
    return keep_lp ? (res,log_probs) : res
end
function assign_taxonomy(seqs::Vector,refs::Vector,taxa::Array; k = 8, n_bootstrap = 100,keep_lp = false,lp=false)
    assignments,log_probs = naieve_bayes(seqs,refs,taxa,k,n_bootstrap,lp)
    res = classification_result(fill("",length(seqs)),seqs,assignments)
    return keep_lp ? (res,log_probs) : res
end
function assign_taxonomy(seq,id,refs,taxa; k = 8, n_bootstrap = 100,keep_lp = false,lp=false)
    assignments,log_probs = naieve_bayes([seq],refs,taxa,k,n_bootstrap,lp)
    res = classification_result(id,seq,assignments)
    return keep_lp ? (res,log_probs) : res
end
function assign_taxonomy(seq,refs,taxa; k = 8, n_bootstrap = 100,keep_lp = false,lp=false)
    assignments,log_probs = naieve_bayes([seq],refs,taxa,k,n_bootstrap,lp)
    res = classification_result("",seq,assignments)
    return keep_lp ? (res,log_probs) : res
end
