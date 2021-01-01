using LegendrePolynomials
using Documenter

DocMeta.setdocmeta!(LegendrePolynomials, :DocTestSetup, :(using LegendrePolynomials); recursive=true)

makedocs(;
    modules=[LegendrePolynomials],
    authors="Jishnu Bhattacharya",
    repo="https://github.com/jishnub/LegendrePolynomials.jl/blob/{commit}{path}#L{line}",
    sitename="LegendrePolynomials",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://jishnub.github.io/LegendrePolynomials.jl",
        assets=String[],
    ),
    pages=[
        "Introduction" => "index.md",
        "Derivatives" => "derivatives.md",
    ],
)

deploydocs(;
    repo="github.com/jishnub/LegendrePolynomials.jl",
)
