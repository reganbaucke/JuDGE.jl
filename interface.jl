workspace();

include("tree.jl")
include("judge.jl")

using JuMP 
using Gurobi

investcost = [180 50 60 40 60 10 10];
investvol = 4;
u_0 = 3;


volume = [ 6   2   1   1  1
    8   2   2   2  1
    8   1   1   1  3
    4   4   3   1  2
    1   3   1   1  2
    7   3   1   1  1
    2   5   2   1  2];

cost = [ 60   20   10   15   10
    8   10   20   20  10
    8   10   15   10  30
    40   40   35   10  20
    15   35   15   15  20
    70   30   15   15  10
    25   50   25   15  20];

mutable struct Knapsack
    itemreward::Array{Float64,1}
    volume::Array{Float64,1}
    p::Float64
    investcost::Float64
end


function P(n::Node)
    list = Node[]
    list = getparents(n)
    push!(list,n)
end

# build a tree
mytree = buildtree(2,2)

getnode(mytree,[1]).data = Knapsack(cost[1,:],volume[1,:],1,investcost[1]);
getnode(mytree,[1,1]).data = Knapsack(cost[2,:],volume[2,:],0.5,investcost[2]);
getnode(mytree,[1,2]).data = Knapsack(cost[3,:],volume[3,:],0.5,investcost[3]);
getnode(mytree,[1,1,1]).data = Knapsack(cost[4,:],volume[4,:],0.25,investcost[4]);
getnode(mytree,[1,1,2]).data = Knapsack(cost[5,:],volume[5,:],0.25,investcost[5]);
getnode(mytree,[1,2,1]).data = Knapsack(cost[6,:],volume[6,:],0.25,investcost[6]);
getnode(mytree,[1,2,2]).data = Knapsack(cost[7,:],volume[7,:],0.25,investcost[7]);


nodes = mytree.nodes
items = 1:5

m = Model(solver=GurobiSolver())

@variable(m, z[nodes], category=:Bin)
capacityvariable(z)

@variable(m, y[nodes,items], category=:Bin)

@objective(m,Min, sum(n.data.p*(n.data.investcost*z[n] - sum(n.data.itemreward[i]*y[n,i] for i in items) ) for n in nodes))

# TODO
# generate the rest of the constraints and tag them appropriatly
@constraint(m, capcon[n in nodes], sum(n.data.volume[i]*y[n,i] for i in items) <= u_0 + sum(investvol*z[h] for h in P(n)))
capacityconstraint(capcon)

@constraint(m,hist[n in nodes],sum(z[h] for h in P(n)) <= 1)
onlyoneinvestment(hist)

solve(m)
println(m.objVal)
println(getvalue(m[:z]))

master  = constructmaster(m)

