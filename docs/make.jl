using NIfTI
using Documenter
using DocumenterVitepress

makedocs(;
    modules=[NIfTI],
    sitename="NIfTI.jl",
    format=DocumenterVitepress.MarkdownVitepress(
        repo = "github.com/JuliaNeuroscience/NIfTI.jl", # this must be the full URL!
        devbranch = "master",
        devurl="dev";
    ),
    pages=[
        "Get Started" => "get_started.md",
        "API" => "api.md"
    ],
)

DocumenterVitepress.deploydocs(;
    repo="github.com/JuliaNeuroscience/NIfTI.jl",
    devbranch = "master",
    push_preview = true,
)
