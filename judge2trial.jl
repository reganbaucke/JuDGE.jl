push!(LOAD_PATH, ".")

using JuDGETree
using JuMP
# using Gurobi
using JuDGE2

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
    bigvol::Float64
    smallvol::Array{Float64,1}
    biginvestcost::Float64
    smallinvestcost::Array{Float64,1}
    volume::Array{Float64,1}
end


function P(n::Node)
    list = Node[]
    list = getparents(n)
    push!(list,n)
end

# build a tree
mytree = buildtree(3,2)

for i = 1:length(mytree.nodes)
    mytree.nodes[i].data = Knapsack(cost[i,:],10,[2,3],90,[40,50],volume[i,:])
end

# mytree[1].data = Knapsack(cost[1,:],10,[2,3],investcost[1]);
# mytree[1,1].data = Knapsack(cost[2,:],volume[2,:],0.5,investcost[2],1:5);
# mytree[1,2].data = Knapsack(cost[3,:],volume[3,:],0.5,investcost[3],1:5);
# mytree[1,1,1].data = Knapsack(cost[4,:],volume[4,:],0.25,investcost[4],1:5);
# mytree[1,1,2].data = Knapsack(cost[5,:],volume[5,:],0.25,investcost[5],1:5);
# mytree[1,2,1].data = Knapsack(cost[6,:],volume[6,:],0.25,investcost[6],1:5);
# mytree[1,2,2].data = Knapsack(cost[7,:],volume[7,:],0.25,investcost[7],1:5);

####
# Set up empty JuDGE Model based off the tree
####
hello = JuDGEModel(mytree)


JuDGEexpansions!(hello) do sp
    # sets
    expa = 1:2
    items = 1:5

    @expansion(sp,small[1:2])
    @expansion(sp,big)
end

JuDGEsubproblem!(hello) do sp, n
    big = sp.objDict[:big]
    small = sp.objDict[:small]

    items = 1:5
    expa = 1:2
    # set up the sub problem variables
    @variable(sp, y[items], category=:Bin)

    #set up objective: note that expansions aren't costed here
    @objective(sp, Min, -n.p * sum(n.data.itemreward[i]*y[i] for i in items))

    # set up the constraints
    @constraint(sp, sum(n.data.volume[i]*y[i] for i in items) <= u_0 + sum(small[e]*n.data.smallvol[e] for e in expa)  + big*n.data.bigvol)
end

# JuDGEexpansioncosts!(hello) do master
    #set up objective: note that expansions aren't costed here
    # @objective(master, Min, -n.p * sum(n.data.itemreward[i]*y[i] for i in items))
# end

tech=[:hydro,:gas]
items = 1:5
expa = 1:2

JuDGEbuild!(hello)
# function JuDGEsubproblems!(f,jmodel::JuDGEModel)
# end
m = Model()
h = @variable(m)
# g = @variable(m,[tech,items,expa])

g = Dict()

for i in items
    g[i] = @variable(m,[Symbol[:hydro,:gas],1:2])
end




# @constraint(m,[t in tech])

some = Dict()
for i in items
    some[i] =  :( @constraint(m,[t in tech, e in expa], 0 <= sum(g[h][t,e] for h in 1:$i )) )
end

ex = some[5]





# ex = Expr(:vect)
# for (i,set) in enumerate(g.indexsets)
#     push!(ex.args,Expr(:call,:in, Symbol(string(i)),set))
# end

# ex2 = Expr(:ref)
# push!(ex2.args,:g)
# for (i,set) in enumerate(g.indexsets)
#     push!(ex2.args, Symbol(string(i)))
# end


# ex2 = :(g[t,e,i])



# ex3 = :(@constraint(m,$(ex),0 <= $(ex2)))
# raft = eval(ex3)
