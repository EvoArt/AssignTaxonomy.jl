using AssignTaxonomy
using Test

@testset "AssignTaxonomy.jl" begin
    ref_fasta = "rdp_train_set_16.fa"
    seq_fasta = "5sp_16S.fasta"
    res = ["Achromobacter", "Achromobacter", "Achromobacter", "Achromobacter", "Ochrobactrum", "Ochrobactrum", 
            "Ochrobactrum", "Ochrobactrum", "Pseudomonas", "Pseudomonas", "Pseudomonas", "Pseudomonas", "Pseudomonas", "Pseudomonas", "Pseudomonas", "Pseudomonas", "Variovorax", "Variovorax", "Variovorax", "Stenotrophomonas", "Stenotrophomonas", "Stenotrophomonas"]
    my_result = assign_taxonomy(seq_fasta,ref_fasta)
    @test all(my_result.Genus .== res)
end
