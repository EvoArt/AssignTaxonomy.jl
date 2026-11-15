using AssignTaxonomy
using Test
const A = AssignTaxonomy
isdefined(A, :setprogress!) && A.setprogress!(false)

randseq(rng, n) = String(rand(rng, ['A', 'C', 'G', 'T'], n))

@testset "AssignTaxonomy.jl" begin
    @testset "R's random number generator" begin
        # set.seed(1); runif(3) and set.seed(42); runif(2) in R
        @test A.runif(A.RMersenneTwister(1), 3) ≈ [0.2655087, 0.3721239, 0.5728534] atol = 2e-7
        @test A.runif(A.RMersenneTwister(42), 2) ≈ [0.9148060, 0.9370754] atol = 2e-7
    end

    @testset "reference parsing" begin
        tax, td = A._parse_lineages([" K;P;C;O;F;G; ", "K;P;C;O;F;", "K;P2;C;O;F;G;"])
        @test td == 6
        @test tax[2] == "K;P;C;O;F;" * A.UNSPECIFIED * ";"
        unite = ["SH$i|AB$i|SH$i.08FU|refs|k__Fungi;p__Asco;c__unidentified;o__O;f__F;g__G;s__G_sp" for i in 1:10]
        tax, td = A._parse_lineages(unite)
        @test tax[1] == "k__Fungi;p__Asco;" * A.UNSPECIFIED * ";o__O;f__F;g__G;" * A.UNSPECIFIED
        @test_throws ErrorException A._parse_lineages(["S001 Genus species"])
    end

    @testset "classifier on a toy reference" begin
        rng = A.Random.MersenneTwister(3)
        base = [randseq(rng, 300) for _ in 1:3]
        lineage = ["Bacteria;Phy1;Cla1;Ord1;Fam1;GenA;", "Bacteria;Phy1;Cla1;Ord1;Fam2;GenB;", "Bacteria;Phy2;Cla2;Ord2;Fam3;GenC;"]
        mutate(s) = (c = collect(s); for p in rand(rng, 1:length(c), 3); c[p] = rand(rng, "ACGT"); end; String(c))
        refs = [mutate(base[g]) for g in 1:3 for _ in 1:5]
        ids = [lineage[g] for g in 1:3 for _ in 1:5]
        model = train_taxonomy(refs, ids)
        q = [mutate(base[2]), A.rc(mutate(base[3])), base[1][1:40]]
        t = assignTaxonomy(q, model; seed=1, outputBootstraps=true)
        @test t.Genus[1] == "GenB"
        @test all(ismissing, t.tax[3, :]) && all(==(0), t.boot[3, :])   # < 50 nt
        @test t.boot[1, 1] == 100
        t2 = assignTaxonomy(q, model; seed=1, tryRC=true)
        @test t2.Genus[2] == "GenC"
        @test assignTaxonomy(q, model; seed=7, outputBootstraps=true).boot ==
              assignTaxonomy(q, model; seed=7, outputBootstraps=true, multithread=false).boot
    end

    @testset "species by exact matching" begin
        hits = Int32[1, 2, 3]
        @test A._maphits(hits, ["coli", "coli", "coli"], 1) == "coli"
        @test ismissing(A._maphits(hits, ["coli", "fergusonii", "coli"], 1))
        @test A._maphits(hits, ["fergusonii", "coli", "coli"], 2) == "coli/fergusonii"
        @test A._maphits(hits, ["Escherichia", "Shigella", "Shigella"], 1) == "Escherichia/Shigella"
        @test ismissing(A._maphits(Int32[], ["x"], 1))
        @test A._maphits(hits, Union{Missing,String}["a", missing, "a"], 2) == "a"
        @test ismissing(A._maphits(hits, Union{Missing,String}["a", missing, "a"], 1))
        # R 4.5 (ICU root collation) sorted these as:
        rsorted = split("_x -x ,x ;x :x !x ?x 'x \"x [x] {x} @x *x &x #x %x `x ^x <x =x >x |x 1a a-b a.b a(b a/b " *
                        "a+b a~b a\$b a1 a10 a2 ab ab aB Ab AB durans ENA07 faecalis sp sP Sp", ' ')
        rsorted = vcat(rsorted, ["sp 1"], split("sp_1 sp-1 sp. x X x_y x9 xa xA xy", ' '))
        @test A._rsort!(String.(reverse(rsorted))) == rsorted
        @test A._matchgenera("Escherichia/Shigella", "Shigella")
        @test A._matchgenera("Clostridium_sensu_stricto_1", "Clostridium")
        @test !A._matchgenera("Clostridiales", "Clostridium")
        @test !A._matchgenera(missing, "X")

        rng = A.Random.MersenneTwister(5)
        r1 = randseq(rng, 400); r2 = randseq(rng, 400)
        mktempdir() do dir
            ref = joinpath(dir, "sp.fa")
            open(ref, "w") do io
                println(io, ">S1 GenA alpha\n", r1)
                println(io, ">S2 GenA beta\n", r1[1:200] * "N" * r1[202:end])
                println(io, ">S3 GenB gamma\n", r2)
            end
            q = [r1[50:250], r1[250:350], A.rc(r2[10:200]), r2[1:100]]
            # q1 spans S2's N (so only S1 holds it); q2 is in S1 and S2; q3 only reverse-complemented
            s = assignSpecies(q, ref)
            @test isequal(s.Genus, ["GenA", "GenA", missing, "GenB"])
            @test isequal(s.Species, ["alpha", missing, missing, "gamma"])
            s = assignSpecies(q, ref; allowMultiple=true, tryRC=true)
            @test isequal(s.Genus, ["GenA", "GenA", "GenB", "GenB"])
            @test isequal(s.Species, ["alpha", "alpha/beta", "gamma", "gamma"])
            tt = TaxonomyTable(q, ["Kingdom", "Genus"], Union{Missing,String}["B" "GenA"; "B" "GenX"; "B" "GenB"; "B" missing],
                               zeros(Int, 4, 0))
            a = addSpecies(tt, ref; tryRC=true)
            @test a.levels == ["Kingdom", "Genus", "Species"]
            @test isequal(a.Species, ["alpha", missing, "gamma", missing])
            @test_throws ErrorException assignSpecies(["ACGTN" * r1[1:60]], ref)
        end
    end

    @testset "assign_taxonomy (original API)" begin
        ref_fasta = "rdp_train_set_16.fa"
        seq_fasta = "5sp_16S.fasta"
        res = ["Achromobacter", "Achromobacter", "Achromobacter", "Achromobacter", "Ochrobactrum", "Ochrobactrum",
               "Ochrobactrum", "Ochrobactrum", "Pseudomonas", "Pseudomonas", "Pseudomonas", "Pseudomonas", "Pseudomonas",
               "Pseudomonas", "Pseudomonas", "Pseudomonas", "Variovorax", "Variovorax", "Variovorax", "Stenotrophomonas",
               "Stenotrophomonas", "Stenotrophomonas"]
        my_result, model = assign_taxonomy(seq_fasta, ref_fasta; keep_lp=true)
        @test all(my_result.Genus .== res)
        @test all(0 .<= my_result.Confidence .<= 1)
        @test model isa TaxonomyModel
        again = assign_taxonomy(seq_fasta, ref_fasta; lp=model)
        @test again.Genus == my_result.Genus
        t = assignTaxonomy(seq_fasta, model; seed=1, minBoot=0)
        @test t.Genus == res
    end
end
