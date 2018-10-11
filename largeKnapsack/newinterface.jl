push!(LOAD_PATH, "..")

using JuDGETree
using JuMP
using Gurobi
using JuDGE2

mutable struct Knapsack
    itemreward::Array{Float64,1}
    volume::Array{Float64,1}
    investcost::Array{Float64,1}
end

(mytree,investvol,initialcap) = deserialize(open("mediumtree.sl"))

investments =  1:2
items = 1:10

m = JuDGE2.JuDGEModel(mytree)

####
# Set up the duals required for each node
####
JuDGEexpansions!(m) do sp
    @expansion(sp,extrabags[1:2])
end

####
# Set up the subproblems for each node
####
JuDGEsubproblem!(m) do sp,n
    # bring expansions into scope for subproblems
    extrabags = sp.objDict[:extrabags]

    items = 1:10

    # set up the sub problem variables
    @variable(sp, y[items], category=:Bin)

    #set up objective: note that expansions aren't costed here
    @objective(sp, Min, -n.p * sum(n.data.itemreward[i]*y[i] for i in items))

    # set up the constraints
    @constraint(sp, sum(n.data.volume[i]*y[i] for i in items) <= initialcap + sum(investvol[o]*extrabags[o] for o in 1:2))
end

####
# Set up the expansion costs
####
JuDGEexpansioncosts!(m) do master,expansion,n
    # bring expansions into scope
    extrabags = expansion[:extrabags]
    return @expression(master,n.p*(sum(n.data.investcost[o]*extrabags[o] for o in 1:2)))
end

JuDGEsolve!(m,10)
# JuDGEpsolve!(m,10,5)

solve(m.master)
println(m.master.objVal)

