using LiebLinigerBA
using Documenter

DocMeta.setdocmeta!(LiebLinigerBA, :DocTestSetup, :(using LiebLinigerBA); recursive=true)

makedocs(;
    modules=[LiebLinigerBA],
    authors="Wei Tang <tangwei@smail.nju.edu.cn> and contributors",
    repo="https://github.com/tangwei94/LiebLinigerBA.jl/blob/{commit}{path}#{line}",
    sitename="LiebLinigerBA.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://tangwei94.github.io/LiebLinigerBA.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/tangwei94/LiebLinigerBA.jl",
)
