push!(LOAD_PATH, "../src/")
push!(LOAD_PATH, "./src/")

using JuDGETree
using JuMP
using Gurobi
using JuDGE

mutable struct Knapsack
    itemreward::Array{Float64,1}
    volume::Array{Float64,1}
    investcost::Array{Float64,1}
end

(mytree,investvol,initialcap) = deserialize(open("./example/mediumtree.sl"))
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
JuDGEsolve!(m) do time, iterations, lb,ub

    println(lb," ",ub)

    if time > 10
        return true
    end
    if iterations > 200
        return true
    end
    if (ub - lb) < 0.001
        return true
    end
    return false
end

####
# get out some of the variables
####
node = mytree.root
finished = false
while !finished
    # this is how you access the value of the binary expansions in the master
    var = m.mastervar[node][:extrabags]

    # print value of variable
    println(node)
    print("Build extra bags?: ")
    print(getvalue(var[1]))
    print(" ")
    println(getvalue(var[2]))

    if !(length(node.children) == 0)
        node = node.children[2]
    else
        finished = true
    end
end