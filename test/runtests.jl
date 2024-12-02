using AssignTaxonomy
using Test


using VectorizedKmers
kmer_array = count_kmers(dna"AACCGGTT", 2)

@testset "AssignTaxonomy.jl" begin
    # Write your tests here.
end

ref_fasta = "rdp_train_set_16.fa"
seq_fasta = "5sp_16S.fasta"

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

function get_targets(seq_fasta)
    ids = []
    seqs = []
    open(FASTA.Reader, seq_fasta) do reader
        for record in reader
                push!(seqs,sequence(record))
                push!(ids,identifier(record))
        end
    end
    seqs,ids
end


word_prior(n,N) = (n+0.5)/(N+1)

function word_priors!(priors ::Array,N)
    for i in eachindex(priors)
        priors[i] = word_prior(priors[i],N)
    end
end

conditional_prob(m,P,M) = (m+P)/(M+1)

function count_mers(refs, k = 8)
    n_refs = length(refs) 
    mer_vec= [Array{Float64}(undef,fill(4,k)...) for _ in 1:n_refs]
    tot_kmer_array= zeros(Float64,fill(4,k)...)
    Threads.@threads for i in 1:n_refs   
        GC.gc()   
        count_ref_mers!(refs[i],mer_vec[i],tot_kmer_array,k)
    end
    return tot_kmer_array, mer_vec
end

function count_ref_mers!(ref_seq,kmer_array,tot_kmer_array,k = 8)
    kmer_array .= count_kmers(LongDNA{4}(ref_seq), k).values .> 0
    tot_kmer_array .+= kmer_array
end

function count_seq_mers(ref_seq,k = 8)
    kmer_array = count_kmers(LongDNA{4}(ref_seq), k).values .> 0
    return kmer_array
end

function assign(seq_mask,log_probs) 
    cond_probs =[sum(log_prob[seq_mask]) for log_prob in log_probs]
    return findmax(cond_probs)[2]
end

function bootstrap(kmer_vec,log_probs,assignment, sample_size,n_bootstrap) 
    hits = 0
    seq_inds = eachindex(kmer_vec)[kmer_vec]
    for i in 1:n_bootstrap
        inds = rand(seq_inds,sample_size) 
        if assign(inds,log_probs)  == assignment
            hits +=1
        end
    end
    return hits/n_bootstrap
end


#### Underlying algorithm
function naieve_bayes(seqs::Vector,refs::Vector,k, n_bootstrap,lp=false)
    t = time()
    N = length(refs)
    n = length(seqs)
    assignments = Vector{Int64}(undef,n)
    confs = Vector{Float64}(undef,n)
    if lp == false
        priors, a =count_mers(refs)
        println(time()-t)
        word_priors!(priors,N) 
        println(time()-t)
        Threads.@threads for i in 1:N
            GC.gc()
            a[i] .= log.(conditional_prob.(a[i],priors,1))
        end
        println(time()-t)
        log_probs = a
        println(time()-t)
    else
        log_probs = lp
    end
    println(time()-t)
    Threads.@threads for seq in seqs
        GC.gc()
        kmer_array = count_seq_mers(seq)
        assignment = assign(kmer_array,log_probs)
        println(time()-t)
        assignments[i] =assignment
            sample_size = sum(kmer_array) ÷ k
        confs[i] = bootstrap(vec(kmer_array),log_probs,assignment,sample_size,n_bootstrap)
        println(time()-t)
    end
    return assignments, confs, log_probs

end


function naieve_bayes(seqs::Vector,refs::Vector,taxa ::Array,k, n_bootstrap,lp=false)
    a,c,l = naieve_bayes(seqs,refs,k, n_bootstrap,lp)
    t = taxa[a,:]
    return hcat(t,c),l
end


#### Top level function

### Most typical function call, taking reference and target fasta input
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
    res = hcat(ids,seqs,assignments)
    return keep_lp ? (res,log_probs) : res
end
function assign_taxonomy(seqs::Vector,refs::Vector,taxa::Array; k = 8, n_bootstrap = 100,keep_lp = false,lp=false)
    assignments,log_probs = naieve_bayes(seqs,refs,taxa,k,n_bootstrap,lp)
    res = hcat(fill("",length(seqs)),seqs,assignments)
    return keep_lp ? (res,log_probs) : res
end
function assign_taxonomy(seq,id,refs,taxa; k = 8, n_bootstrap = 100,keep_lp = false,lp=false)
    assignments,log_probs = naieve_bayes([seq],refs,taxa,k,n_bootstrap,lp)
    res = hcat(id,seq,assignments)
    return keep_lp ? (res,log_probs) : res
end
function assign_taxonomy(seq,refs,taxa; k = 8, n_bootstrap = 100,keep_lp = false,lp=false)
    assignments,log_probs = naieve_bayes([seq],refs,taxa,k,n_bootstrap,lp)
    res = hcat("",seq,assignments)
    return keep_lp ? (res,log_probs) : res
end
ass = assign_taxonomy(seq_fasta,ref_fasta)

@benchmark assign_taxonomy(LongDNA{4}(seqs[1]),"baccy",ref_fasta,n_bootstrap = 20)
@profview assign_taxonomy(seq_fasta,ref_fasta,n_bootstrap = 50,lp =L)

@profview assign_taxonomy(seq_fasta,ref_fasta,lp = l)