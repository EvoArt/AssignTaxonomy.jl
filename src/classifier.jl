
#https://doi.org/10.1128/AEM.00062-07
#The Classifier uses a feature space consisting of all possible 8-base
#subsequences (words).
# The psition of the word in a sequence is ignored.
# Only the words occurring in the query contribute to the score.

base_vals = Base.ImmutableDict(DNA_A=>0,DNA_C=>1,DNA_G=>2,DNA_T=>3)

word_prior(n,N) = (n+0.5)/(N+1)
function word_priors!(priors ::Vector,N)
    for i in 1:N
        priors[i] = word_prior(priors[i],N)
    end
end

conditional_prob(m,P,M) = (m+P)/(M+1)

function get_mer_idx(mer)
    idx = 0
    for base in mer
        idx = 4*idx + base_vals[base]
    end
    return idx +1
end


function count_ref_mers!(ref_seq,column_idx,mer_array,mer_vec,k = 8)
    for mer in each(DNAMer{k}, ref_seq)
        current_mer = mer.fw
        current_mer_idx = get_mer_idx(current_mer)
        mer_array[current_mer_idx,column_idx] = 1
        mer_vec[current_mer_idx] +=1
    end
end

function count_seq_mers(seq,k = 8)
    mask = zeros(Bool,4^k)
    for mer in each(DNAMer{k}, seq)
        current_mer = mer.fw
        current_mer_idx = get_mer_idx(current_mer)
        mask[current_mer_idx] = 1
    end
    return mask, (1:4^k)[mask]
end


function count_mers(refs, k = 8)
    n_refs = length(refs) 
    mer_vec= zeros(Float64,4^k)
    mer_array= zeros(Bool,4^k,n_refs)
    for i in 1:n_refs
        count_ref_mers!(refs[i],i,mer_array,mer_vec,k)
    end
    return mer_vec, mer_array
end

function count_mers(refs, genus_vec ::Vector, k = 8)
    n_refs = length(refs) 
    mer_vec= zeros(Float64,4^k)
    genus_array= zeros(Int64,4^k,maximum(genus_vec))
    mer_array= zeros(Bool,4^k,n_refs)
    for i in 1:n_refs
        count_ref_mers!(refs[i],i,mer_array,mer_vec,k)
        genus_array[:,genus_vec[i]] .+= mer_array[:,i]
    end
    return mer_vec, genus_array
end


function assign(seq_mask,log_probs) 
    cond_probs = vec(sum(log_probs[seq_mask,:], dims = 1))
    println(maximum(cond_probs))
    return findmax(cond_probs)[2]
end

function bootstrap(seq_inds,log_probs,assignment, sample_size) 
    hits = 0
    for i in 1:100
        inds = rand(seq_inds,sample_size) 
        cond_probs = vec(sum(log_probs[inds,:], dims = 1))
        if findmax(cond_probs)[2] == assignment
            hits +=1
        end
    end
    return hits/100
end

function naieve_bayes(seqs,refs)
    N = length(refs)
    assignments = []
    confs = []
    priors, a =count_mers(refs)
    word_priors!(priors,N) 
    log_probs = log.(conditional_prob.(a,priors,1)) 
    for seq in seqs
        seq_mask, seq_inds = count_seq_mers(seq)
        sample_size = sum(seq_mask) ÷8
        assignment = assign(seq_mask,log_probs)
        push!(assignments,assignment)
        push!(confs, bootstrap(seq_inds,log_probs,assignment,sample_size))
    end
    return assignments, confs
end

function naieve_bayes(seqs,refs,genus_vec)
    N = length(refs)
    assignments = []
    confs = []
    M = reverse([sum(genus_vec .==i) for i in 1:maximum(genus_vec)]) # reverse gives right answer. but not sure why
    priors, a =count_mers(refs,genus_vec)
    word_priors!(priors,N) 
    log_probs = log.(conditional_prob.(a,priors,M')) #M` gave wrong results. need to figure out right way...
    for seq in seqs
        seq_mask, seq_inds = count_seq_mers(seq)
        sample_size = sum(seq_mask)  ÷8
        assignment = assign(seq_mask,log_probs)
        push!(assignments,assignment)
        push!(confs, bootstrap(seq_inds,log_probs,assignment,sample_size))
    end
    return assignments, confs
end

function get_genera(fasta)
    genera =[]
    open(FASTA.Reader,fasta) do reader
        for record in reader
            push!(genera,split(identifier(record),";")[5])
        end
    end
    unq = unique(genera)
    n = length(unq)
    mapping = Base.ImmutableDict([unq[i] => i for i in 1:n]...)
    reverse_mapping = Base.ImmutableDict([i => unq[i] for i in 1:n]...)
    return [mapping[genus] for genus in genera], reverse_mapping
end

function get_genera(genera ::Vector)
    unq = unique(genera)
    n = length(unq)
    mapping = Base.ImmutableDict([unq[i] => i for i in 1:n]...)
    reverse_mapping = Base.ImmutableDict([i => unq[i] for i in 1:n]...)
    return [mapping[genus] for genus in genera], reverse_mapping
end

function assignTaxonomy(seqs,ref_fasta)
    genera = []
    refs = []
    open(FASTA.Reader, ref_fasta) do reader
        for record in reader
            if length(findall(";",identifier(record))) >4
                push!(refs,sequence(record))
                push!(genera,split(identifier(record),";")[6])
            end
        end
    end
    genus_vec, genus_mapping = get_genera(genera)

    assignments, confs = naieve_bayes(seqs,refs,genus_vec)
    return hcat(assignments,[genus_mapping[a] for a in assignments],confs)
    
end


using FASTX
