using NIfTI
using Documenter

makedocs(;
    modules=[NIfTI],
    repo="https://github.com/JuliaNeuroscience/NIfTI.jl/blob/{commit}{path}#{line}"
    sitename="NIfTI.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaNeuroscience.github.io/NIfTI.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaNeuroscience/NIfTI.jl",
)
