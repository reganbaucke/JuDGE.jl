using Documenter, JuDGE

makedocs(sitename="JuDGE - Julia Decomposition for General Expansion",
        modules = [JuDGE],
        format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"),
        pages = ["Main" => "index.md"])
