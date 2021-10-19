using AssignTaxonomy
using Documenter

DocMeta.setdocmeta!(AssignTaxonomy, :DocTestSetup, :(using AssignTaxonomy); recursive=true)

makedocs(;
    modules=[AssignTaxonomy],
    authors="Arthur Newbury",
    repo="https://github.com/EvoArt/AssignTaxonomy.jl/blob/{commit}{path}#{line}",
    sitename="AssignTaxonomy.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://EvoArt.github.io/AssignTaxonomy.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/EvoArt/AssignTaxonomy.jl",
)
