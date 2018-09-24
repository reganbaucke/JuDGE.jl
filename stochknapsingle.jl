using Gurobi
using JuMP


# describe data
# sets
nodes = 1:7
items = 1:5

# parameters
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

# invest
investcost = [180 50 60 40 60 20 60];
investvolume = [4 4 4 4 4 4 4];

# probabilities of nodes
prob = [ 1.00 0.50 0.50 0.25 0.25 0.25 0.25];

# initial capacity
u_0 = 3;

# initialize predecesor array
P = Array{Array{Int,1},1}(7)

P[1] = [1];
P[2] = [1, 2];
P[3] = [1, 3];
P[4] = [1, 2, 4];
P[5] = [1, 2, 5];
P[6] = [1, 3, 6];
P[7] = [1, 3, 7];


m = Model(solver=GurobiSolver(OutputFlag=0))

# @variable(m, 0 <= y[items,nodes] <=1 )
# @variable(m, 0 <= z[nodes] <=1)

@variable(m, y[items,nodes], category=:Bin )
@variable(m, z[nodes],  category=:Bin)

@objective(m,Min, sum(prob[n]*(investcost[n] * z[n] - sum(cost[n,i] * y[i,n] for i in items)) for n in nodes) )

con =[]
for n in nodes
    push!(con,@constraint(m, sum(volume[n,i]*y[i,n] for i in items) <= u_0 + sum(investvolume[sn]*z[sn] for sn in P[n])))
    @constraint(m, sum(z[sn] for sn in P[n]) <= 1)
end

solve(m)
println(m.objVal)
println(getvalue(m[:z]))


# println(getvalue(m[:z]))
# println(getvalue(m[:y]))
# getdual(con[3])