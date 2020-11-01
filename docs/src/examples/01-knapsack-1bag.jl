# # Knapsack Problem - Single Expansion
# This first example demonstrates how a basic JuDGE model can be set up, with
# a 7-node tree, and a single expansion.
using JuMP, JuDGE, GLPK

# First we will define our tree with a depth of 2 and a degree of 2.
mytree = narytree(2,2)

# ## Data

# Now we specify data for each node of the tree. Here we define the data using
# functions, but you could use dictionaries, as shown in subsequent examples.
function invest_cost(node)
   if node == get_node(mytree,[1])
	  180.0
   elseif node == get_node(mytree,[1,1])
	  50.0
   elseif node == get_node(mytree,[1,2])
	  60.0
   elseif node == get_node(mytree,[1,1,1])
	  40.0
   elseif node == get_node(mytree,[1,1,2])
	  60.0
   elseif node == get_node(mytree,[1,2,1])
	  10.0
   elseif node == get_node(mytree,[1,2,2])
	  10.0
   end
end

function item_volume(node)
   if node == get_node(mytree,[1])
	  [6, 2, 1, 1, 1]
   elseif node == get_node(mytree,[1,1])
	  [8, 2, 2, 2, 1]
   elseif node == get_node(mytree,[1,2])
	  [8, 1, 1, 1, 3]
   elseif node == get_node(mytree,[1,1,1])
	  [4, 4, 3, 1, 2]
   elseif node == get_node(mytree,[1,1,2])
	  [1, 3, 1, 1, 2]
   elseif node == get_node(mytree,[1,2,1])
	  [7, 3, 1, 1, 1]
   elseif node == get_node(mytree,[1,2,2])
	  [2, 5, 2, 1, 2]
   end
end

function item_reward(node)
   if node == get_node(mytree,[1])
	  [60, 20, 10, 15, 10]
   elseif node == get_node(mytree,[1,1])
	  [8, 10, 20, 20, 10]
   elseif node == get_node(mytree,[1,2])
	  [8, 10, 15, 10, 30]
   elseif node == get_node(mytree,[1,1,1])
	  [40, 40, 35, 10, 20]
   elseif node == get_node(mytree,[1,1,2])
	  [15, 35, 15, 15, 20]
   elseif node == get_node(mytree,[1,2,1])
	  [70, 30, 15, 15, 10]
   elseif node == get_node(mytree,[1,2,2])
	  [25, 50, 25, 15, 20]
   end
end

# ## Subproblem definitions

# The JuDGE subproblems are defined through a function that takes a node as its
# single argument, returning a JuMP model. We use GLPK as our optimizer.
JuDGE_SP_Solver = optimizer_with_attributes(GLPK.Optimizer, "msg_lev" => 0, "mip_gap" => 0.0)
function sub_problems(node)
   model = Model(JuDGE_SP_Solver)
   @expansion(model, bag)
   @expansioncosts(model, bag*invest_cost(node))
   @variable(model, y[1:5], Bin)
   @constraint(model, BagExtension, sum(y[i]*item_volume(node)[i] for i in 1:5) <= 3 + 4 * bag)
   @objective(model, Min, sum(-item_reward(node)[i] * y[i] for i in 1:5))
   model
end

# Note that there are a few JuDGE-specific macros used to define a JuDGE subproblem

# This defines our single expansion variable
@expansion(model, bag)

# This defines the cost of buying this bag at particular node
@expansioncosts(model, bag*invest_cost(node))

# JuDGE subproblem objective is defined using the standard @objective macro:
@objective(model, Min, sum(-item_reward(node)[i] * y[i] for i in 1:5))


# ## Defining and solving the JuDGE model
# The JuDGEModel is defined based on a tree, probability distribution,
# sub_problems, and an optimizer using an interior point method.
JuDGE_MP_Solver = optimizer_with_attributes((method=GLPK.INTERIOR) -> GLPK.Optimizer(),
							"msg_lev" => 0, "mip_gap" => 0.0)
judy = JuDGEModel(mytree, ConditionallyUniformProbabilities, sub_problems, JuDGE_MP_Solver)

# Solve the model
JuDGE.solve(judy)

# ## Displaying the output
# Print the objective, and optimal expansions
println("Objective: "*string(judy.bounds.UB))
JuDGE.print_expansions(judy,onlynonzero=false)

# Re-solve the subproblems and print the objective
println("Re-solved Objective: " * string(resolve_subproblems(judy)))

# ## Deterministic equivalent
# Here we set up and solve the deterministic equivalent, this uses the same structure as a `JuDGEModel`.
# We set up our GLPK solver to show its progress.
JuDGE_DE_Solver = optimizer_with_attributes(GLPK.Optimizer, "msg_lev" => 2, "mip_gap" => 0.0)
deteq = DetEqModel(mytree, ConditionallyUniformProbabilities, sub_problems, JuDGE_DE_Solver)
JuDGE.solve(deteq)
println("Deterministic Equivalent Objective: " * string(objective_value(deteq.problem)))
