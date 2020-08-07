# Tutorials

## Tutorial 1: A basic JuDGE model

### Problem description
For our tutorial, we will consider the following optimization problem: Our goal
is minimize the cost of a stochastic sequence of knapsack problems. We represent
the stochastic process with a discrete scenario tree. At each node in the
scenario tree, we solve a knapsack problem. However at any point in the tree, we
have the ability to expand the capacity of our knapsack at a certain cost. Once
the capacity of our bag has been expanded, we are able use the extra volume for
future knapsack problems from this node of the tree forward.

In this optimization problem, we are trading off against the cost of expanding
our knapsack, versus the ability to fit more into our knapsack. Deciding when to
perform the knapsack expansion is the difficult part of this
optimization problem.


### Solving our problem using JuDGE

First things first, we should bring into our workspace JuDGE, JuMP, and Gurobi using
```@example tutorial
using JuDGE
using JuMP
using GLPK
```
Using these libraries, we will now model and solve our optimization problem.


The lifecycle of a `JuDGEModel` is the following:

1. The definition of a `JuDGETree`;
2. defining the subproblems of the `JuDGEModel`;
3. building the `JuDGEModel`;
3. solving the `JuDGEModel`.

The user's job is both Steps 1 and 2, while JuDGE will automatically perform
Steps 3 and 4.

A `JuDGETree` can be built in many different ways. A `JuDGETree` simply consists
of the root node of the tree, and a list of all the nodes in the tree. This is
defined as a nested set of subtrees, with the final nodes being leaf nodes. Each
subtree simply defines its children, and there are functions that facilitate the
calculation of its parent and the probability of arriving at the node, and the
data that correspondes to the node, can be referenced through dictionaries.

For now, we will build a tree of depth 3, where each node has 2 children with
uniform probabilities using `buildtree`:
```@example tutorial
mytree = narytree(2,2)
```
`mytree` is a tree which contains 7 nodes, with depth 2, and degree 2.
(A depth of 0, gives only a single leaf node.) We can visualise the tree using
```@example tutorial
JuDGE.print_tree(mytree)
```
### Problem data
Let us now associate our tree with the problem data.

For our instance of the problem, we will use the following data: for each node,
we will have a knapsack problem with three items to choose from, each with
different rewards, and different volumes. The structure of this data is arbitrary;
JuDGE just needs to be able to access the relevant data, based on the node being
processed (dictionaries or functions are recommended).
```@example tutorial
invest_cost = Dict( zip( collect(mytree,order=:breadth), [15, 8, 8, 4, 4, 4, 4]) )


item_volume = Dict( zip( collect(mytree,order=:breadth), [ [4, 3, 3, 1, 2],
                                                           [5, 3, 4, 2, 1],
                                                           [5, 4, 2, 7, 2],
                                                           [5, 4, 1, 8, 2],
                                                           [3, 1, 5, 6, 3],
                                                           [2, 5, 8, 4, 6],
                                                           [7, 5, 4, 2, 3] ]) )

item_reward = Dict( zip( collect(mytree,order=:breadth), [ [32, 9, 9, 4, 8],
                                                           [30, 12, 40, 10, 9],
                                                           [20, 28, 12, 42, 12],
                                                           [40, 28, 9, 24, 10],
                                                           [15, 7, 20, 48, 12],
                                                           [10, 30, 54, 32, 30],
                                                           [32, 25, 24, 14, 24] ]) )
```
We can print the tree again and display a parameter for each node:
```@example tutorial
JuDGE.print_tree(mytree,item_reward)
```
We also define some other parameters that apply to all the nodes.
```@example tutorial
num_items=5
num_invest=6
initial_volume = 6
invest_volume = [2,2,2,3,3,3]
```
### Subproblems
We need to now define the subproblems. These are JuMP models with some JuDGE-
specific features. For our knapsack problem:
```@example tutorial
JuDGE_SP_Solver = optimizer_with_attributes(GLPK.Optimizer, "msg_lev" => 0, "mip_gap" => 0.0)
function sub_problems(node)
   sp = Model(JuDGE_SP_Solver)
   @expansion(sp, invest[1:num_invest])
   @capitalcosts(sp, sum(invest[i]*invest_volume[i] for i=1:num_invest)*invest_cost[node])
   @variable(sp, y[1:num_items], Bin)
   @constraint(sp, BagExtension, sum(y[i]*item_volume[node][i] for i in 1:num_items) <=
        initial_volume + sum(invest_volume[i] * invest[i] for i in 1:num_invest))
   @sp_objective(sp, sum(-item_reward[node][i] * y[i] for i in 1:num_items))
   return sp
end
```
The three elements of this that make it a JuDGE subproblem are:

`@expansion(model, bag)` This defines the expansion variables, and supports
standard JuMP vectorized variable declaration. These will be binary.

`@capitalcosts` This declares an expression for the costs of investment; this
must be linear (an AffExpr).

The overall optimization problem at each node problem is a classical knapsack
problem. We have specified the initial volume of the knapsack is `initial_volume`, and
the each expansion increases the volume of the knapsack by the corresponding value in
`invest_volume`.

### Solving the JuDGE Model
Now with our tree built and the problem data referenced, we can initialize the
`JuDGEModel` based on our tree, subproblems, and solver.
```@example tutorial
JuDGE_MP_Solver = optimizer_with_attributes((method=GLPK.INTERIOR) -> GLPK.Optimizer(),
							"msg_lev" => 0, "mip_gap" => 0.0)
judy = JuDGEModel(mytree, ConditionallyUniformProbabilities, sub_problems, JuDGE_MP_Solver)
```
`ConditionallyUniformProbabilities` simply applies a uniform conditional probability
distribution for child nodes. Any function mapping nodes to absolute probabilities
can be used here.

At this point, we have constructed a valid `JuDGEModel`.
We can now solve our model by making a call to `JuDGE.solve`:
```@example tutorial
JuDGE.solve(judy)
```
There are a number of optional stopping critieria that can be set here:
    `abstol`, `reltol`, `rlx_abstol`, `rlx_reltol`, `duration`, `iter`.

Currently, we recommend using JuDGE with Gurobi as the subproblem and master problem
solvers. Any solvers can be specified, but the master problem must return duals, and
a barrier method is recommended to computational efficiency. The subproblems can be
solved with any method, but currently need to be solved to optimality (bound gap of 0).

We can view the optimal solution to our problem by calling
```@example tutorial
println("Objective: "*string(judy.bounds.UB))
JuDGE.print_expansions(judy)
```
Finally, if we want to recover the optimal solutions for the nodes, we must fix the
investments and resolve each subproblem, after which we can write the solution to
a CSV file.
```@example tutorial
println("Re-solved Objective: " * string(resolve_subproblems(judy)))
```
```
JuDGE.write_solution_to_file(judy,joinpath(@__DIR__,"knapsack_solution.csv"))
```
## Tutorial 2: Formatting output
Since the expansion variables are all binary, the `print_expansions` function doesn't
directly convey how much capacity is built. In order to output more information about the capacity
it's possible to write a custom function that is supplied as an optional argument `format`.
First, let us simply provide the capacity associated with each decision.
```@example tutorial
function format_output(s::Symbol,value)
   if s==:invest
      output=Dict{Int64,Float64}()
      for i in 1:num_invest
         output[i]=invest_volume[i]*value[i]
      end
	  return output
   end
   return nothing
end

JuDGE.print_expansions(judy, format=format_output)
```
We can also aggregate the capacities of all the expansions.
```@example tutorial
function format_output(s::Symbol,value)
   if s==:invest
	  return sum(invest_volume[i]*value[i] for i in 1:num_invest)
   end
   return nothing
end

JuDGE.print_expansions(judy, format=format_output)
```
## Tutorial 3: Ongoing costs
Depnding on the capacity planning application that JuDGE is being applied to there may be
ongoing upkeep / maintenance costs for the expansions. This is modelled within JuDGE by
using the `@ongoingcosts` macro to specify the cost of the expansions being available
at each node in the scenario tree. (That is, the corresponding expansion variable has
been set to 1 in the master problem.) If there are costs that can be avoided, by not making
use of capacity that is granted, then this should be modelled within the `@sp_objective`.

```@example tutorial
function sub_problems_ongoing(node)
   sp = Model(JuDGE_SP_Solver)
   @expansion(sp, invest[1:num_invest])
	@capitalcosts(sp, sum(invest[i]*invest_volume[i] for i=1:num_invest)*invest_cost[node])
	@ongoingcosts(sp, sum(invest[i]*invest_volume[i] for i=1:num_invest)*2)
   @variable(sp, y[1:num_items], Bin)
   @constraint(sp, BagExtension, sum(y[i]*item_volume[node][i] for i in 1:num_items) <=
        initial_volume + sum(invest_volume[i] * invest[i] for i in 1:num_invest))
   @sp_objective(sp, sum(-item_reward[node][i] * y[i] for i in 1:num_items))
   return sp
end

judy = JuDGEModel(mytree, ConditionallyUniformProbabilities, sub_problems_ongoing,
	JuDGE_MP_Solver)

JuDGE.solve(judy)

println("Objective: "*string(judy.bounds.UB))
JuDGE.print_expansions(judy, format=format_output)
```
## Tutorial 4: Deterministic equivalent
JuDGE can use the tree, and subproblems to automatically construct the deterministic equivalent of the
stochastic capacity expansion problem. This is created by defining a `DetEqModel` with the same arguments
as a `JuDGEModel`, as follows:
```@example tutorial
JuDGE_DE_Solver = optimizer_with_attributes(GLPK.Optimizer, "msg_lev" => 2, "mip_gap" => 0)
deteq = DetEqModel(mytree, ConditionallyUniformProbabilities, sub_problems, JuDGE_DE_Solver)
JuDGE.solve(deteq)
```
The solution can be printed to the REPL or a CSV in the same way as a `JUDGEModel`.
```@example tutorial
println("Deterministic Equivalent Objective: " * string(objective_value(deteq.problem)))
JuDGE.print_expansions(deteq, format=format_output)
```
```
JuDGE.write_solution_to_file(judy,joinpath(@__DIR__,"knapsack_solution.csv"))
```
## Tutorial 5: Lag and duration
Depending on the particular expansion problem being modelled, there may be some delay (lag)
between when the expansion decision is made, and when the capacity becomes available.
This can be modelled in JuDGE be specifying a lag when defining the expansion variables in the
subproblems, the `@expansion` marco is overloaded to allow an additional argument which is an
integer specifying the lag, e.g. for a lag of one you would define the expansion variables as follows:
```
@expansion(sp, invest[1:num_invest],1)
```
The duration of an expansion is set to ∞ by default. However, if an expansion is temporary or otherwise has a limited
lifespan, we can include an additional argument when defining the expansion variable. For example if the lag is 0 and the
duration is 2, we can define it as follows:
```
@expansion(sp, invest[1:num_invest],0,2)
```
We will now redefine our subproblems and re-solve our model. (Note that the `invest_cost` has been divided by 2.)
```@example tutorial
function sub_problems_lag(node)
   sp = Model(JuDGE_SP_Solver)
   @expansion(sp, invest[1:num_invest],1)
   @capitalcosts(sp, sum(invest[i]*invest_volume[i] for i=1:num_invest)*invest_cost[node]/2)
   @variable(sp, y[1:num_items], Bin)
   @constraint(sp, BagExtension, sum(y[i]*item_volume[node][i] for i in 1:num_items) <=
        initial_volume + sum(invest_volume[i] * invest[i] for i in 1:num_invest))
   @sp_objective(sp, sum(-item_reward[node][i] * y[i] for i in 1:num_items))
   return sp
end

judy = JuDGEModel(mytree, ConditionallyUniformProbabilities, sub_problems_lag, JuDGE_MP_Solver)
JuDGE.solve(judy)

println("Objective: "*string(judy.bounds.UB))
println("Lower Bound: "*string(judy.bounds.LB))
JuDGE.print_expansions(judy, format=format_output)
```
We see that although this solution has ostensibly converged, the best integer is greater than the lower bound.
The * in the final line of the output means that the integer solution has been found using a MIP solve for the
generated columns. In order to find a provably optimal solution we must use branch and price.

## Tutorial 6: Branch and Price
JuDGE implements a branch and price algorithm for problems which are not naturally integer. It can be run with default settings
as follows:
```@example tutorial
judy = JuDGEModel(mytree, ConditionallyUniformProbabilities, sub_problems_lag, JuDGE_MP_Solver)
best = JuDGE.branch_and_price(judy,search=:lowestLB,branch_method=JuDGE.constraint_branch)

println("Objective: "*string(best.bounds.UB))
println("Lower Bound: "*string(best.bounds.LB))
JuDGE.print_expansions(best, format=format_output)
```
We now see that we have found a better solution, and proved it is optimal.

There are several options for the search: `:lowestLB` always chooses to branch on the node with the lowest lower bound;
`:depth_first_dive` performs a depth-first search of the branch and bound tree, once it find a node with an integer relaxation
it keeps adjacent nodes within the tree; `:depth_first_resurface` performs a depth-first search of the branch and bound tree,
but once it find a node with an integer relaxation it returns to the root node and explores the other branch; `:breadth_first`
performs a breadth-first search of the tree.

There are two built-in branching methods: `JuDGE.constraint_branch` and `JuDGE.variable_branch`. It is also possible
to write custom branching methods; see the API for more details.

## Tutorial 7: Risk aversion
JuDGE implements risk aversion using the risk measure CVaR over the accumulated profits up to each of leaf nodes in the scenario tree.
The objective function minimized is a convex combination of expectation and CVaR, with the parameter λ=1 meaning at all the weight is
placed on CVaR. In our implementation CVaR represents the expected cost of the 100β% worst outputs. In order to implement CVaR, we supply
the optional argument `CVaR=(λ,β)` when we construct the `JuDGEModel`.
```@example tutorial
judy = JuDGEModel(mytree, ConditionallyUniformProbabilities, sub_problems, JuDGE_MP_Solver,
	CVaR=(0.5,0.1))
best = JuDGE.branch_and_price(judy,search=:lowestLB,branch_method=JuDGE.constraint_branch)

println("Objective: "*string(best.bounds.UB))
println("Lower Bound: "*string(best.bounds.LB))
JuDGE.print_expansions(best, format=format_output)
```

## Tutorial 8: Shutdown variables
JuDGE supports shutdown decisions. These are variables that remove capacity when they are activated. The `@capitalcosts` of these
decisions may be negative, reflecting some salvage value; moreover, the `@ongoingcosts` may also be negative, reflecting avoided
maintenance costs. Given that this is a shutdown variable, it is important to remember that the capacity being removed should be accounted
for elsewhere within the subproblem.

```@example tutorial
function sub_problems_shutdown(node)
   model = Model(JuDGE_SP_Solver)
   @shutdown(model, bag)
   @capitalcosts(model, -bag*invest_cost[node])
   @variable(model, y[1:num_items], Bin)
   @constraint(model, BagExtension, sum(y[i]*item_volume[node][i] for i in 1:num_items) <=
		7 - bag)
   @sp_objective(model, sum(-item_reward[node][i] * y[i] for i in 1:num_items))
   return model
end

judy = JuDGEModel(mytree, ConditionallyUniformProbabilities, sub_problems_shutdown,
	JuDGE_MP_Solver)
JuDGE.solve(judy)

println("Objective: "*string(objective_value(judy.master_problem)))
JuDGE.print_expansions(judy)
```

## Tutorial 9: Side-constraints
JuDGE supports side-constraints being added to the master problem. These can be
constraints across expansion variables at a single node, or can be constraints on
variables corresponding to different nodes. The `JuDGE.history` function can be
useful if applying logical constraints on expansion variables. In order to apply
a budget constraint at each node we can define the following function:
```@example tutorial
function budget(model,tree)
   for node in collect(tree)
      @constraint(model,sum(invest_cost[node]*invest_volume[i]*invest[node][i]
		  for i in 1:num_invest)<=40)
   end
end
```
We now can define a `JuDGEModel` with these side-constraints, and solve it using branch and price.
```@example tutorial
judy = JuDGEModel(mytree, ConditionallyUniformProbabilities, sub_problems, JuDGE_MP_Solver,
	sideconstraints=budget)
judy = JuDGE.branch_and_price(judy,search=:lowestLB,branch_method=JuDGE.constraint_branch)
JuDGE.print_expansions(judy, format=format_output)
```
