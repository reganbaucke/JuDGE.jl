workspace();
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
    25   50   25   15  20]

# invest
investcost = [180 50 60 40 60 10 10];
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


# set up master problem
master = Model(solver=GurobiSolver(OutputFlag=0,Method=2))

@variable(master, 0 <= x[nodes] <=1)
@objective(master,Min, sum(prob[n]*investcost[n]*x[n] for n in nodes))


for n in nodes
    @constraint(master, sum(x[sn] for sn in P[n]) <= 1)
end

# set up empty constraints
pi = Array{JuMP.ConstraintRef{JuMP.Model,JuMP.GenericRangeConstraint{JuMP.GenericAffExpr{Float64,JuMP.Variable}}},1}(length(nodes))
for n in nodes
    pi[n] = @constraint(master, 0 <= sum(x[h] for h in P[n]))
end

# for n in nodes
#     @constraint(master, sum(x[h] for h in P[n]) <=1)
# end

mu = Array{Any,1}(length(nodes))
for n in nodes
    mu[n] =  @constraint(master, 0 == 1)
end

# μ = similar(collect(nodes))
μ = zeros(length(nodes))

# π_n = Array{Float64,1}(length(nodes))
π = zeros(length(nodes))

# set up node problems
sub = []
for n in nodes
    # subprob= Model(solver=GurobiSolver())
    subprob = Model(solver=GurobiSolver(OutputFlag=0))
    @variable(subprob, z, category=:Bin)
    @variable(subprob, y[items], category=:Bin)
    @objective(subprob, Min, prob[n]*sum(-cost[n,i]*y[i] for i in items) - π[n]*z - μ[n])
    @constraint(subprob, sum(volume[n,i]*y[i] for i in items) <= u_0 + investvolume[n]*z)
    push!(sub,subprob)
end

solve(master)

# populate the duals
π = getdual(pi);
μ = getdual(mu);

#with these duals, solve the sub problems
for n in nodes
    @objective(sub[n], Min, prob[n]*sum(-cost[n,i]*sub[n][:y][i] for i in items) - π[n]*sub[n][:z] - μ[n])
    solve(sub[n])
end

#for each node somehow create a the column coeffecients
# first get the objective coefficents
objcoef = Array{Float64,1}(length(nodes));
for n in nodes
    objcoef[n] = prob[n]*sum(-cost[n,i]*getvalue(sub[n][:y][i]) for i in items)
end

colcoef = Array{Float64,1}(length(nodes))
for n in nodes
    colcoef[n] = getvalue(sub[n][:z])
end

# add each column into the master problem
for n in nodes
    # @variable(master, 0 <= w <= 1, objective = objcoef[n], inconstraints = [pi[n];mu[n]], coefficients = [colcoef;1])
    # @variable(master, 0 <= w <= 1, objective = objcoef[n], inconstraints = pi[n], coefficients = colcoef[n])
    Variable(master,0.0, 1.0,:Cont, objcoef[n], [pi[n];mu[n]], [colcoef[n];1],"w_$n^1")
end


currentObj = -99999
for j = 1:20
    solve(master)

    if abs(master.objVal- currentObj) <= 0.0001
        println("----------------------")
        println("took $j iterations")
        println("objective function value: $(master.objVal)")
        println(getvalue(master[:x]))
        break;
    else
        currentObj = master.objVal
    end

    # populate the duals
    π = getdual(pi)
    μ = getdual(mu);

    #with these duals, solve the sub problems
    for n in nodes
        @objective(sub[n], Min, prob[n]*sum(-cost[n,i]*sub[n][:y][i] for i in items) - π[n]*sub[n][:z] - μ[n])
        # solve(sub[n])
    end

    for n in nodes
        # @objective(sub[n], Min, prob[n]*sum(-cost[n,i]*sub[n][:y][i] for i in items) - π[n]*sub[n][:z] - μ[n])
        solve(sub[n])
    end

    #for each node somehow create a the column coeffecients
    # first get the objective coefficents
    objcoef = Array{Float64,1}(length(nodes));
    for n in nodes
        objcoef[n] = prob[n]*sum(-cost[n,i]*getvalue(sub[n][:y][i]) for i in items)
    end

    colcoef = Array{Float64,1}(length(nodes))
    for n in nodes
        colcoef[n] = getvalue(sub[n][:z])
    end

    # add each column into the master problem
    for n in nodes
        # @variable(master, 0 <= w <= 1, objective = objcoef[n], inconstraints = [pi[n];mu[n]], coefficients = [colcoef;1])
        # @variable(master, 0 <= w <= 1, objective = objcoef[n], inconstraints = pi[n], coefficients = colcoef[n])
        Variable(master,0.0, 1.0,:Cont, objcoef[n], [pi[n];mu[n]], [colcoef[n];1],"w_$n^$(j+1)")
    end
end

