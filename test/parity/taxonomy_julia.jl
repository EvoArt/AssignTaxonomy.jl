# Julia half of the taxonomy parity check (see taxonomy_r.R). Exits 1 on a
# difference that is not explained by an exact tie between genera (R breaks
# those with a non-reproducible random draw).
#   julia --project=. test/parity/taxonomy_julia.jl <outdir> <train_set.fa.gz> <species.fa.gz>
using AssignTaxonomy
const A = AssignTaxonomy
A.setprogress!(false)

rdir, ref, spref = ARGS
f(name) = joinpath(rdir, name)
seqs = readlines(f("queries.txt"))
ex_ref, ex_sp = readlines(f("extdata.txt"))
println("R collation locale: ", readchomp(f("locale.txt")))

nfail = 0
function check(name, ok, detail="")
    global nfail += !ok
    println(ok ? "PASS " : "FAIL ", name, ok || isempty(detail) ? "" : "  -- " * detail)
end
function readtab(path)
    L = readlines(path)
    rows = [split(l, '\t'; keepempty=true) for l in L[2:end]]
    return split(L[1], '\t'), [[x == "NA" ? missing : String(x) for x in r] for r in rows]
end
cells(t::TaxonomyTable, i) = [t.tax[i, l] for l in 1:size(t.tax, 2)]

runs = Dict("default" => (ref, false, 100, 50), "rc" => (ref, true, 100, 50),
            "boot80" => (ref, false, 7, 80), "example" => (ex_ref, false, 1, 50))
models = Dict{String,TaxonomyModel}()
for nm in ["default", "rc", "boot80", "example"]
    r, tryRC, seed, minBoot = runs[nm]
    t0 = time()
    model = get!(() -> train_taxonomy(r), models, r)
    t1 = time()
    t = assignTaxonomy(seqs, model; minBoot, tryRC, seed, outputBootstraps=true)
    nm == "default" && println("  Julia ($(Threads.nthreads()) threads): train $(round(t1 - t0; digits=1)) s, ",
                               "classify $(round(time() - t1; digits=1)) s; R (1 thread, train + classify): ",
                               readchomp(f("time_default.txt")), " s")
    _, _, ties = A._classify(model, seqs; tryRC, seed)
    hdr, rtax = readtab(f("tax_$nm.tsv"))
    _, rboot = readtab(f("boot_$nm.tsv"))
    check("taxLevels $nm", t.levels == hdr)
    taxbad = [i for i in eachindex(seqs) if !isequal(cells(t, i), rtax[i])]
    bootbad = [i for i in eachindex(seqs) if t.boot[i, :] != parse.(Int, rboot[i])]
    unexplained = [i for i in union(taxbad, bootbad) if !ties[i]]
    detail = "tax differs for $(length(taxbad)), boot for $(length(bootbad)); not explained by ties: $unexplained"
    bd = [sum(abs.(t.boot[i, :] .- parse.(Int, rboot[i]))) for i in bootbad]
    println("  $nm: $(count(ties)) of $(length(seqs)) sequences met an exact tie; ",
            "$(length(taxbad)) differ in taxonomy and $(length(bootbad)) in bootstraps (all tied), ",
            "bootstrap |difference| summed over ranks: median $(isempty(bd) ? 0 : sort(bd)[cld(length(bd), 2)]), max $(maximum(bd; init=0))")
    for i in first(unexplained, 3)
        println("    seq $i: J ", cells(t, i), " ", t.boot[i, :], "\n           R ", rtax[i], " ", rboot[i])
    end
    check("assignTaxonomy $nm (tax and bootstraps)", isempty(unexplained), detail)
end

acgt = readlines(f("acgt.txt")) .== "1"
q = seqs[acgt]
for (file, sp, kw) in (("sp_single.tsv", spref, (;)), ("sp_multi_rc.tsv", spref, (; allowMultiple=true, tryRC=true)),
                       ("sp_multi3.tsv", spref, (; allowMultiple=3)), ("sp_example.tsv", ex_sp, (; allowMultiple=true)))
    t = assignSpecies(q, sp; kw...)
    hdr, r = readtab(f(file))
    bad = [i for i in eachindex(q) if !isequal(cells(t, i), r[i])]
    for i in first(bad, 3)
        println("    seq $i: J ", cells(t, i), "  R ", r[i])
    end
    check("assignSpecies $file", t.levels == hdr && isempty(bad), "$(length(bad)) differ")
end

let (hdr, rt) = readtab(f("tax_default.tsv"))        # R's table as input (isolated)
    tax = Matrix{Union{Missing,String}}(permutedims(reduce(hcat, rt[acgt])))
    sub = TaxonomyTable(q, String.(hdr), tax, zeros(Int, count(acgt), 0))
    for (file, kw) in (("addsp.tsv", (;)), ("addsp_multi_rc.tsv", (; allowMultiple=true, tryRC=true)))
        t = addSpecies(sub, spref; kw...)
        hdr, r = readtab(f(file))
        bad = [i for i in eachindex(q) if !isequal(cells(t, i), r[i])]
        check("addSpecies $file", t.levels == hdr && isempty(bad), "$(length(bad)) differ")
    end
end

println(nfail == 0 ? "ALL PASS" : "$nfail FAILED")
exit(nfail == 0 ? 0 : 1)
