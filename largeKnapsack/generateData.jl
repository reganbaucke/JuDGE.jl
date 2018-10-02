workspace()
push!(LOAD_PATH, "..")

using JuDGETree
using JuMP
using Gurobi
using JuDGE

# how many investments?
numinvest = 4;

# number of items to pick from in the knapsack?
numitems = 20

# size of tree?
degree = 4
depth = 8

totalnodes = Int64((degree^depth-1)/(degree-1))

investcost = zeros(totalnodes,numinvest)
for i = 1:totalnodes
    investcost[i,:] = rand(numinvest)*2 + [5.5,6.5,7.5,8.5]
end
 
investvol = [40,45,50,70]
initialcap = 80

itemvolume = zeros(totalnodes,numitems)
for i = 1:totalnodes
    itemvolume[i,:] = ((rand(numitems)-0.5)*2)*2 + collect(linspace(4,22,numitems))
end

itemcost = zeros(totalnodes,numitems)
for i = 1:totalnodes
    itemcost[i,:] = ((rand(numitems)-0.5)*2)*0.5 + collect(linspace(0.5,1,numitems))
end

mutable struct Knapsack
    itemreward::Array{Float64,1}
    volume::Array{Float64,1}
    p::Float64
    investcost::Array{Float64,1}
end

mytree = buildtree(depth-1,4)

for i in 1:totalnodes
    mytree.nodes[i].data =  Knapsack(itemcost[i,:], itemvolume[i,:], 0.0, investcost[i,:] )
end

wholetree!(mytree) do n
    if n == mytree.root
        n.data.p = 1
    else
        n.data.p = n.parent.data.p/degree
    end
end

open("bigtree.sl", "w") do file
    serialize(file,(mytree,investvol,initialcap,investcost))
end

