using Documenter, ECON622_BLP
makedocs(;
         modules=[ECON622_BLP],
         format=Documenter.HTML(),
         pages=[
           "Home" => "index.md",
           "Function Reference" => "functions.md"
         ],
         repo=Remotes.GitHub("UBCECON","ECON622_BLP.jl"),
         sitename="ECON622_BLP.jl",
         #doctest=false,
         authors="Paul Schrimpf <schrimpf@mail.ubc.ca>"
         #assets=String[],
         #plugins=[bib]
)


deploydocs(repo="github.com/UBCECON/ECON622_BLP.jl.git", branch="gh-pages")
