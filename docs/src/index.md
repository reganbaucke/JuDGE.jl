```@meta
CurrentModule = JuDGE
DocTestSetup = quote
    using JuDGE
end
```
![JuDGE](assets/judge-small.png)

# Getting started


## Requirements

JuDGE requires Julia-1.3+, JuMP and appropriate optimiser(s). For academics,
Gurobi / CPLEX provide free academic licenses, otherwise, you can use CBC/Clp or
GLPK.

## Installation

JuDGE is installed by the `Pkg` utility provided by Julia. In the Julia REPL,
simply make the following function call.

    ] add "https://github.com/reganbaucke/JuDGE.jl"

Then, in your Julia script, use

    using JuDGE
to import the functions from the JuDGE module into the current namespace.
