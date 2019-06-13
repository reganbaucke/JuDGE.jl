using JuMP
using Gurobi
using JuDGE

mutable struct Knapsack
    itemreward::Array{Float64,1}
    volume::Array{Float64,1}
    investcost::Array{Float64,1}
end

(mytree,investvoltmp,initialcap) = deserialize(open("C:/Users/adow031/.julia/v0.6/JuDGE/example/mediumtree.sl"))

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

function printNodeExpansionsDet(node::Node,onlynonzero::Bool)
    # this is how you access the value of the binary expansions in the master

    # print value of variable (if non-zero)
    output=false
    if onlynonzero
        if getvalue(z[node,1])>0 || getvalue(z[node,2])>0
            output=true
        end
    else
        output=true
    end
    if output
        print("Expansion decisions: ")
        println(node)
        print("(1) ")
        println(getvalue(z[node,1]))
        print("(2) ")
        println(getvalue(z[node,2]))
    end

    if !(length(node.children) == 0)
        for i in 1:length(node.children)
            printNodeExpansionsDet(node.children[i],onlynonzero)
        end
    end
end

printNodeExpansionsDet(mytree.root,true)
