using ElectricFields
using Documenter

isdefined(Main, :NOPLOTS) && NOPLOTS || include("plots.jl")

makedocs(;
    modules=[ElectricFields],
    authors="Stefanos Carlström <stefanos.carlstrom@gmail.com> and contributors",
    sitename="ElectricFields.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://www.tipota.org/ElectricFields.jl",
        assets=String[],
        mathengine = MathJax(Dict(
            :TeX => Dict(
                :equationNumbers => Dict(:autoNumber => "AMS"),
                :Macros => Dict(
                    :defd => "≝",
                    :abs => ["\\left|#1\\right|", 1],
                    :vec => ["\\mathbf{#1}", 1],
                    :mat => ["\\mathsf{#1}", 1],
                    :conj => ["#1^*", 1],
                    :im => "\\mathrm{i}",
                    :diff => ["\\mathrm{d}#1\\,", 1],
                    :bmat => ["\\begin{bmatrix}#1\\end{bmatrix}", 1]
                )
            )
        ))
    ),
    pages=[
        "Home" => "index.md",
        "Interface" => "interface.md",
        "Field types" => "field_types.md",
        "Envelopes" => "envelopes.md",
        "Carriers" => "carriers.md",
        "Field properties" => "properties.md",
        "Dispersion" => "dispersion.md",
        "Tabulated fields" => "tabulated_fields.md",
        "Reference" => "reference.md",
    ],
    doctest=false,
    checkdocs=:exports
)

deploydocs(;
    repo="github.com/jagot/ElectricFields.jl",
    push_preview = true,
)
