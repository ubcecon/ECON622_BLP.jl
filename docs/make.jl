using Documenter, ECON622_BLP
makedocs(;
         modules=[ECON622_BLP],
         format=Documenter.HTML(),
         pages=[
           "Home" => "index.md",
           "Function Reference" => "functions.md"
         ],
         repo="https://github.com/UBCECON/ECON622_BLP.jl/blob/{commit}{path}#L{line}",
         sitename="ECON622_BLP.jl",
         #doctest=false,
         authors="Paul Schrimpf <schrimpf@mail.ubc.ca>"
         #assets=String[],
)


deploydocs(repo="github.com/UBCECON/ECON622_BLP.jl.git")
