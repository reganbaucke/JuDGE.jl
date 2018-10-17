push!(LOAD_PATH, "../src/")

using JuDGETree
using JuMP
using Gurobi
using JuDGE

mutable struct Knapsack
    itemreward::Array{Float64,1}
    volume::Array{Float64,1}
    investcost::Array{Float64,1}
end

(mytree,investvoltmp,initialcap) = deserialize(open("./example/mediumtree.sl"))

investments =  1:length(investvoltmp)
items = 1:10

investvol = investvoltmp


det = Model(solver=GurobiSolver())

@variable(det, y[mytree.nodes,items], category=:Bin )
@variable(det, z[mytree.nodes,investments],  category=:Bin)

@objective(det,Min, sum(n.p*(sum(n.data.investcost[o]/2 * z[n,o] for o in investments) - sum(n.data.itemreward[i] * y[n,i] for i in items)) for n in mytree.nodes))

@constraint(det,[n in mytree.nodes], sum(n.data.volume[i]*y[n,i] for i in items) <= initialcap + sum(sum(investvol[o]*z[h,o] for h in P(n)) for o in investments))
@constraint(det,[n in mytree.nodes, o in investments], sum(z[h,o] for h in P(n)) <= 1)

JuMP.solve(det)

optimalValue = det.objVal