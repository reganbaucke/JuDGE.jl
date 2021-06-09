using Test

"""
A list of files to ignore when running the tests.
"""
const _EXCLUDE_EXAMPLES = [
    "advanced.jl",
    "custom_trees.jl",
    "vrp.jl", # too hard / random
    "cutting_stock.jl", # random
    "inventory.jl", # random
    "knapsacks.jl", # random
    "transportation.jl", # delimitedfiles
]

const _EXAMPLES_DIR = joinpath(dirname(@__DIR__), "examples")

# Setup and initialize the subproblem solvers for the example using GLPK.
include(joinpath(_EXAMPLES_DIR, "solvers", "setup_glpk.jl"))

for file in readdir(_EXAMPLES_DIR)
    if !endswith(file, ".jl") || file in _EXCLUDE_EXAMPLES
        continue
    end
    @testset "$(file)" begin
        include(joinpath(_EXAMPLES_DIR, file))
    end
end
