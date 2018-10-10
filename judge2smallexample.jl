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
    volume::Array{Float64,1}
    investcost::Float64
    items::UnitRange{Int64}
end

# build a tree
# mytree = JudgeTree.buildtree(2,2)
mytree = buildtree(3,2)

mytree[1].data = Knapsack(cost[1,:],volume[1,:],investcost[1],1:5);
mytree[1,1].data = Knapsack(cost[2,:],volume[2,:],investcost[2],1:5);
mytree[1,2].data = Knapsack(cost[3,:],volume[3,:],investcost[3],1:5);
mytree[1,1,1].data = Knapsack(cost[4,:],volume[4,:],investcost[4],1:5);
mytree[1,1,2].data = Knapsack(cost[5,:],volume[5,:],investcost[5],1:5);
mytree[1,2,1].data = Knapsack(cost[6,:],volume[6,:],investcost[6],1:5);
mytree[1,2,2].data = Knapsack(cost[7,:],volume[7,:],investcost[7],1:5);


####
# Set up empty JuDGE Model based off the tree
####
hello = JuDGEModel(mytree)

####
# Set up the sets
####
items = 1:5
tech = [:gas, :hydro]

####
# Set up the expansion variables
####
JuDGEexpansions!(hello) do sp
    @expansion(sp,bag)
end

####
# Set up the subproblems for each node
####
JuDGEsubproblem!(hello) do sp, n
    # bring expansions into scope for subproblems
    bag = sp.objDict[:bag]

    items = 1:5

    # set up the sub problem variables
    @variable(sp, y[items], category=:Bin)

    #set up objective: note that expansions aren't costed here
    @objective(sp, Min, -n.p * sum(n.data.itemreward[i]*y[i] for i in items))

    # set up the constraints
    @constraint(sp, sum(n.data.volume[i]*y[i] for i in items) <= u_0 + bag*investvol)
end

####
# Set up the objectives. We do this seperatley from the model as we recall this function several times
####
JuDGEexpansioncosts!(hello) do master,expansion,n
    # bring expansions into scope
    bag = expansion[:bag]

    return @expression(master,n.p*(n.data.investcost*bag))
end

JuDGEbuild!(hello)

JuDGEsolve!(hello,10)

# println(hello.master)
# solve(hello.master)
# println(hello.master.objVal)