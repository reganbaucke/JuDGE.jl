using Documenter, JuDGE, Literate

for dir in ["examples"]
    for file in sort(readdir(joinpath(@__DIR__, "src", dir)))
        !endswith(file, ".jl") && continue
        filename = joinpath(@__DIR__, "src", dir, file)
        Literate.markdown(filename, dirname(filename); documenter=true)
    end
end

const EXAMPLES = Any[
    "examples/$(file)" for file  in filter(
        f -> endswith(f, ".md"),
        sort(readdir(joinpath(@__DIR__, "src", "examples")))
    )
]

makedocs(sitename="JuDGE - Julia Decomposition for General Expansion",
        modules = [JuDGE],
        clean = true,
        format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"),
        pages = ["JuDGE" => "index.md",
                 "JuDGE API" => "judge.md",
                 "AbstractTree API" => "tree.md",
                 "Examples" => EXAMPLES])
