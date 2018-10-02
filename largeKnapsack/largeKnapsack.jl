push!(LOAD_PATH, "..")

using JuDGETree
using JuMP
using Gurobi
using JuDGE

mutable struct Knapsack
    itemreward::Array{Float64,1}
    volume::Array{Float64,1}
    p::Float64
    investcost::Array{Float64,1}
end

thistree = deserialize(open("bigtree.sl"))