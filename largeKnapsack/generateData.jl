push!(LOAD_PATH, "..")

using JuDGETree
using JuMP
using Gurobi
using JuDGE

# how many investments?
numinvest = 2;

# number of items to pick from in the knapsack?
numitems = 20

# size of tree?
degree = 3
depth = 6

totalnodes = Int64((degree^depth-1)/(degree-1))

investcost = zeros(totalnodes,numinvest)
for i = 1:totalnodes
    # investcost[i,:] = rand(numinvest)*2 + [5.5,6.5,7.5,8.5]
    investcost[i,:] = (rand(numinvest)*2  + 2*[2.0,3.5])*(1-((i-1)/(totalnodes*1.2)))
end


 
# investvol = [40,45,50,70]
investvol = [40,50]
initialcap = 80

itemvolume = zeros(totalnodes,numitems)
for i = 1:totalnodes
    # itemvolume[i,:] = ((rand(numitems)-0.5)*2)*2 + collect(linspace(4,22,numitems))
    itemvolume[i,:] = ((rand(numitems)-0.5)*2)*2 + collect(linspace(4,22,numitems))
end

itemcost = zeros(totalnodes,numitems)
for i = 1:totalnodes
    itemcost[i,:] = ((rand(numitems)-0.5)*2)*0.5 + collect(linspace(0.5,1,numitems))
    itemcost[i,:] = ((rand(numitems)-0.5)*2)*0.5 + collect(linspace(0.5,1,numitems))
end

mutable struct Knapsack
    itemreward::Array{Float64,1}
    volume::Array{Float64,1}
    investcost::Array{Float64,1}
end

mytree = buildtree(depth,degree)

for i in 1:totalnodes
    mytree.nodes[i].data =  Knapsack(itemcost[i,:], itemvolume[i,:], investcost[i,:] )
end

open("mediumtree.sl", "w") do file
    serialize(file,(mytree,investvol,initialcap))
end

