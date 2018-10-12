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
m = JuDGE2.JuDGEModel(mytree)

####
# Set up the expansions variables in the problem
####
JuDGEexpansions!(m) do sp
    @expansion(sp,extrabags[1:2])
end

####
# Set up the subproblems for each node
####
JuDGEsubproblem!(m) do sp,n,expansion
    # bring expansions into scope for subproblems, this is just useful for write the JuMP constraints etc.
    extrabags = expansion[:extrabags]

    items = 1:10

    # set up the sub problem variables
    @variable(sp, y[items], category=:Bin)

    #set up objective: note that expansions aren't costed here
    @objective(sp, Min, -n.p * sum(n.data.itemreward[i]*y[i] for i in items))

    # set up the constraints
    @constraint(sp, sum(n.data.volume[i]*y[i] for i in items) <= initialcap + sum(investvol[o]*extrabags[o] for o in 1:2))
end

####
# Set up the expansion costs. This function must return an expression.
####
JuDGEexpansioncosts!(m) do master,n,expansion
    # bring expansions into scope
    extrabags = expansion[:extrabags]

    return @expression(master,(sum(n.data.investcost[o]*extrabags[o] for o in 1:2)))
end

####
# Solve problem with a particular convergence test.
####
JuDGEsolve!(m) do time, iterations, lb,ub

    println(lb," ",ub)

    if time > 10
        return true
    end
    if iterations > 10
        return true
    end
    return false
end