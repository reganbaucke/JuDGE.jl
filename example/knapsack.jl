using JuDGETree
using JuMP
#using Cbc
using Gurobi
using JuDGE

mutable struct Knapsack
    itemreward::Array{Float64,1}
    volume::Array{Float64,1}
    investcost::Array{Float64,1}
end

(mytree,investvol,initialcap) = deserialize(open(joinpath(@__DIR__,"mediumtree.sl")))
m = JuDGEModel(mytree)

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
    @objective(sp, Min, sum(-n.data.itemreward[i]*y[i] for i in items))

    # set up the constraints
    @constraint(sp, sum(n.data.volume[i]*y[i] for i in items) <= initialcap + sum(investvol[o]*extrabags[o] for o in 1:2))
end

####
# Set up the expansion costs. This function must return an expression.
####
JuDGEexpansioncosts!(m) do master,n,expansion
    # bring expansions into scope
    extrabags = expansion[:extrabags]

    return @expression(master,0.5*(sum(n.data.investcost[o]*extrabags[o] for o in 1:2)))
end

####
# Solve problem with a particular convergence test.
####
JuDGEsolve!(m,GurobiSolver(OutputFlag=0)) do time, iterations, lb,ub
    println(lb," ",ub)

    if time > 10
        return true
    end
    if iterations > 200
        return true
    end
    if (ub - lb) < 0.00001
        return true
    end
    return false
end

####
# Solve problem as a deterministic equivalent.
####
JuDGEsolvedeteq!(m,GurobiSolver(OutputFlag=1))
#JuDGEsolvedeteq!(m,CbcSolver())

####
# find the optimal expansion variables
####
printExpansions(m)
printDetEqExpansions(m)
