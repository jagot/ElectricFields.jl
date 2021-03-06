using ElectricFields
using Documenter

makedocs(;
    modules=[ElectricFields],
    authors="Stefanos Carlström <stefanos.carlstrom@gmail.com> and contributors",
    repo="https://github.com/jagot/ElectricFields.jl/blob/{commit}{path}#L{line}",
    sitename="ElectricFields.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://jagot.github.io/ElectricFields.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jagot/ElectricFields.jl",
)
