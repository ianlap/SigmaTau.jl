using Documenter
using SigmaTau, SigmaTauBase, SigmaTauStability, SigmaTauEnsemble

DocMeta.setdocmeta!(SigmaTauBase,      :DocTestSetup, :(using SigmaTau); recursive=true)
DocMeta.setdocmeta!(SigmaTauStability, :DocTestSetup, :(using SigmaTau); recursive=true)
DocMeta.setdocmeta!(SigmaTauEnsemble,  :DocTestSetup, :(using SigmaTau); recursive=true)

makedocs(
    sitename = "SigmaTau.jl",
    modules  = [SigmaTau, SigmaTauBase, SigmaTauStability, SigmaTauEnsemble],
    format   = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical  = "https://ianlap.github.io/SigmaTau.jl",
        mathengine = Documenter.KaTeX(),
    ),
    pages = [
        "Home" => "index.md",
    ],
    doctest  = true,
    warnonly = [:missing_docs, :cross_references],
)

deploydocs(
    repo         = "github.com/ianlap/SigmaTau.jl.git",
    push_preview = true,
    devbranch    = "main",
)
