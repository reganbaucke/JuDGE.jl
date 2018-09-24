include("tree.jl")
using JuMP 
using Gurobi

investcost = [180 50 60 40 60 20 60];
investvol = 4;


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

mutable struct JuDGEModel
    tree::Tree
    modeldesc::JuMP.Model
    master::JuMP.Model
    node::JuMP.Model
end

function P(n::Node)
    list = Node[]
    list = getparents(n)
    push!(list,n)
end

# build a tree
mytree = buildtree(2,2)

getnode(mytree,[1]).data = Knapsack(cost[1,:],volume[1,:],1,180);
getnode(mytree,[1,1]).data = Knapsack(cost[2,:],volume[2,:],0.5,50);
getnode(mytree,[1,2]).data = Knapsack(cost[3,:],volume[3,:],0.5,60);
getnode(mytree,[1,1,1]).data = Knapsack(cost[4,:],volume[4,:],0.25,40);
getnode(mytree,[1,1,2]).data = Knapsack(cost[5,:],volume[5,:],0.25,60);
getnode(mytree,[1,2,1]).data = Knapsack(cost[6,:],volume[6,:],0.25,20);
getnode(mytree,[1,2,2]).data = Knapsack(cost[7,:],volume[7,:],0.25,60);



nodes = mytree.nodes
items = 1:5

m = Model(solver=GurobiSolver())

@variable(m, z[nodes], category=:Bin)
@variable(m, y[nodes,items], category=:Bin)

m.ext[:z] = Int[];
# tag these variables as the capacity variables
for n in nodes
    push!(m.ext[:z],m.objDict[:z][n].col)
end

@objective(m,Min, sum(n.data.p*(n.data.investcost*z[n] - sum(n.data.itemreward[i]*y[n,i] for i in items) ) for n in nodes))


# TODO
# generate the rest of the constraints and tag them appropriatly
con =[]
n = nodes[1]
@constraint(m, sum(volume[n,i]*y[i,n] for i in items) <= u_0 + sum(investvolume[sn]*z[sn] for sn in P(n)))

# for n in nodes
    # push!(con,@constraint(m, sum(volume[n,i]*y[i,n] for i in items) <= u_0 + sum(investvolume*z[sn] for sn in P(n)))
    # @constraint(m, sum(z[sn] for sn in P(n) <= 1)
# end