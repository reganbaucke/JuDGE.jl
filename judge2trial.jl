push!(LOAD_PATH, ".")

using JuDGETree
using JuMP
using Gurobi
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


JuDGEexpansions!(hello) do sp,n
    # set up the expansion variables
    #put the expansion symbols into the sp.ext part of the model
    expa = 1:2
    items = 1:5

    @expansion(sp,small[e in 1:2])
    @expansion(sp,big)

    # set up the sub problem variables
    @variable(sp, y[items], category=:Bin)

    # set up the constraints
    @constraint(sp, sum(n.data.volume[i]*y[i] for i in items) <= u_0 + sum(small[e]*n.data.smallvol[e] for e in expa) + big*n.data.bigvol)

    #set up objective: note that expansions aren't costed here
    @objective(sp, Min, -n.p * sum(n.data.itemreward[i]*y[i] for i in items))
end

# function JuDGEsubproblems!(f,jmodel::JuDGEModel)

# end

function thisadd(a,b,c)
    a+b+c
end

ex = Expr(:call,:thisadd,1,2,3)

ex = Expr(:call,:=,a,3)