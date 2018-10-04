push!(LOAD_PATH, "..")

using JuDGETree
using JuMP
using Gurobi
using JuDGE

mutable struct Knapsack
    itemreward::Array{Float64,1}
    volume::Array{Float64,1}
    p::Float64
    investcost::Array{Float64,1}
end

function P(n::Node)
    list = Node[]
    list = getparents(n)
    push!(list,n)
end

(mytree,investvoltmp,initialcap,investcosttmp) = deserialize(open("bigtree.sl"))

investments =  1:4
items = 1:20

investvol = investvoltmp
investcost = investcosttmp


m = JuDGEModel(mytree)

####
# Set up the duals required for each node
####
JuDGEduals!(m) do n
    @dual(pi,investments)
    @dual(mu)
    return [pi,mu]
end

####
# Set up the subproblems for each node
####
JuDGEsubproblems!(m) do n
    # Set up an empty jump model
    sp = JuMP.Model(solver=GurobiSolver(OutputFlag=0))

    # Set up the variables
    @variable(sp, z[investments], category=:Bin)
    @variable(sp, y[items], category=:Bin)

    # set up the constraints
    @constraint(sp, sum(n.data.volume[i]*y[i] for i in items) <= initialcap + sum(investvol[o]*z[o] for o in investments))

    return sp 
end

####
# Set up the objectives. We do this seperatley from the model as we recall this function several times
####
JuDGEobjective!(m) do n
    # bring into scope to make writing the objective easier
    sp = m.subprob[n]
    dual = m.duals[n]
    data = n.data

    @objective(sp, Min, 
        -data.p * sum(data.itemreward[i]*sp[:y][i] for i in items) - 
        sum(dual[:pi][o]*sp[:z][o] for o in investments) - 
        dual[:mu][])
end

####
# Write out the master problem
####
JuDGEmaster!(m) do
    master = Model(solver=GurobiSolver(OutputFlag=0,Method=2))

    @variable(master, 0 <= x[n in m.tree.nodes, o in investments] <=1)
    @objective(master,Min, sum(n.data.p*sum(n.data.investcost[o]*x[n,o] for o in investments) for n in m.tree.nodes))

    # give these constraints a name so that we can get the duals out of them later
    @constraint(master, pi[n in m.tree.nodes, o in investments] ,0 <= sum(x[h,o] for h in P(n)))
    @constraint(master, mu[n in m.tree.nodes], 0 == 1)

    return master
end

####
# populate the duals from the solution of the master problem
####
JuDGEupdateduals!(m) do n
    # hello.master
    for o in investments
        m.duals[n][:pi][o] = getdual(m.master.objDict[:pi][n,o])
    end
    m.duals[n][:mu][] = getdual(m.master.objDict[:mu][n])
end

####
# tell JuDGE how to build a column from a nodal solution
####
JuDGEbuildcolumn!(m) do n 
    sp = m.subprob[n]
    lb = 0;
    ub = 1;
    obj = n.data.p*sum(-n.data.itemreward[i]*getvalue(sp[:y][i]) for i in items)
    contr = [m.master.objDict[:pi][n,:]; m.master.objDict[:mu][n]]
    colcoef = Array{Float64,1}()
    for o = 1:length(investments)
        push!(colcoef,getvalue(sp[:z][o]))
    end
    push!(colcoef,1.0)
    name = "w"

    return (lb,ub,:Cont ,obj, contr, colcoef, name)
end

JuDGE.solve(m,20)