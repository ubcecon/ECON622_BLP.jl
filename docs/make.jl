using Documenter, ECON622_BLP
makedocs(;
    modules=[ECON622_BLP],
    format=Documenter.HTML(),
    authors="Paul Schrimpf <paul.schrimpf@ubc.ca> and contributors",
    pages=[
        "Home" => "index.md",
        "Function Reference" => "functions.md"
    ],
    repo=Remotes.GitHub("UBCECON","ECON622_BLP.jl"),
    sitename="ECON622_BLP.jl",
    #doctest=false,
    #assets=String[],
    #plugins=[bib]
)

deploydocs(;repo="ubcecon/ECON622_BLP.jl", branch="gh-pages")
