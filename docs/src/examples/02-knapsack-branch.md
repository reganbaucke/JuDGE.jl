```@meta
EditURL = "<unknown>/F0PKY/docs/src/examples/02-knapsack-branch.jl"
```

# Knapsack Problem - Branch and Price
This example demonstrates how the branch and price algorithm can be used to
solve a JuDGE model.

```@example 02-knapsack-branch
using JuMP, JuDGE, GLPK, Random
```

We will randomly generate data for this problem, but to ensure a consistent
result we first set the random seed.

```@example 02-knapsack-branch
Random.seed!(50)
```

We specify that there are 5 possible expansion investments, and 10 items to
select from at each node.

```@example 02-knapsack-branch
numinvest = 5
numitems = 10
```

The tree is defined to have a depth and 4 and a degree of 3.

```@example 02-knapsack-branch
mytree = narytree(4,3)
```

Find the number of nodes in the tree

```@example 02-knapsack-branch
totalnodes = JuDGE.count(mytree)
```

## Data
Set up the investment cost data

```@example 02-knapsack-branch
investcost = zeros(totalnodes,numinvest)
for i = 1:totalnodes
  investcost[i,:] = ([1,1.8,3.5,6.8,13.5])*(1-((i-1)/(totalnodes*1.2)))
end
```

Specify the investment volumes using a binary expansion, with an initial capacity
of 0

```@example 02-knapsack-branch
investvol = [1,2,4,8,16]
initialcap = 0
```

Randomly specify the item volumes and costs

```@example 02-knapsack-branch
itemvolume = zeros(totalnodes,numitems)
for i = 1:totalnodes
  itemvolume[i,:] = ((rand(numitems))*2) + collect(range(4,22,length = numitems))
end

itemcost = zeros(totalnodes,numitems)
for i = 1:totalnodes
  itemcost[i,:] = ((rand(numitems) .- 0.5)*2)*2# + collect(range(0.5,1,length = numitems))
end
```

Create a function that maps the node, and data name to a value

```@example 02-knapsack-branch
nodes = collect(mytree)
function data(node, input)
  input[findall(x -> x == node, nodes)[1], :]
end
```

## Subproblem definitions
Define the subproblems as a function of the node

```@example 02-knapsack-branch
JuDGE_SP_Solver = optimizer_with_attributes(GLPK.Optimizer, "msg_lev" => 0, "mip_gap" => 0.0)
function sub_problems(node)
  model = Model(JuDGE_SP_Solver)
  @expansion(model, bag[1:numinvest])
  @expansioncosts(model, sum(data(node,investcost)[i] * bag[i] for i in  1:numinvest))
  @variable(model, y[1:numitems], Bin)
  @constraint(model, BagExtension ,sum( y[i]*data(node,itemvolume)[i] for i in 1:numitems) <=
  				initialcap + sum(bag[i]*investvol[i] for i in 1:numinvest))
  @sp_objective(model, sum(-data(node,itemcost)[i] * y[i] for i in 1:numitems))
  model
end
```

## Defining and solving the JuDGE model
Set up the JuDGE model.

```@example 02-knapsack-branch
JuDGE_MP_Solver = optimizer_with_attributes((method=GLPK.INTERIOR) -> GLPK.Optimizer(),
							"msg_lev" => 0, "mip_gap" => 0.0)
judy = JuDGEModel(mytree, ConditionallyUniformProbabilities, sub_problems, JuDGE_MP_Solver)
```

Solve the JuDGE model using the branch and price algorithm, using a
constraint branch and branching on the lowest lower bound nodes first

```@example 02-knapsack-branch
best=JuDGE.branch_and_price(judy,rlx_abstol=10^-7,inttol=10^-6,
				branch_method=JuDGE.constraint_branch,search=:lowestLB)
```

## Displaying the output
Print the objective function value

```@example 02-knapsack-branch
println("Objective: "*string(best.bounds.UB))
```

Create a custom function that formats the output for print_expansions(...),
and print the expansions using this function

```@example 02-knapsack-branch
function format_output(s::Symbol,values)
  if s==:bag
	 return sum(values[i]*investvol[i] for i in 1:numinvest)
  end
  return nothing
end

JuDGE.print_expansions(best,format=format_output)
```

## Deterministic equivalent
Set up and solve the deterministic equivalent problem

```@example 02-knapsack-branch
JuDGE_DE_Solver = optimizer_with_attributes(GLPK.Optimizer, "msg_lev" => 2, "mip_gap" => 2)
deteq = DetEqModel(mytree, ConditionallyUniformProbabilities, sub_problems, JuDGE_DE_Solver)
JuDGE.solve(deteq)
println("Deterministic Equivalent Objective: " * string(objective_value(deteq.problem)))
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

