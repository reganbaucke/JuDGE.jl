# # Knapsack Problem - Single Expansion
# This first example demonstrates how a basic JuDGE model can be set up, with
# a 7-node tree, and a single expansion.

using JuMP, JuDGE

# First we will define our tree
mytree = narytree(2,2)

# ### Data

# Now we specify data for each node of the tree.

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

# ### Subproblems

# The JuDGE subproblems are defined through a function that takes a node as its
# single argument, returning a JuMP model.

function sub_problems(node)
   model = Model(JuDGE_SP_Solver)
   @expansion(model, bag)
   @expansioncosts(model, bag*invest_cost(node))
   @variable(model, y[1:5], Bin)
   @constraint(model, BagExtension, sum(y[i]*item_volume(node)[i] for i in 1:5) <= 3 + 4 * bag)
   @sp_objective(model, sum(-item_reward(node)[i] * y[i] for i in 1:5))
   model
end

# Note that there are a few JuDGE-specific macros used to define a JuDGE subproblem

# This defines our single expansion variable
@expansion(model, bag)

# This defines the cost of buying this bag at particular node
@expansioncosts(model, bag*invest_cost(node))

# Instead of using @objective, JuDGE models use:
@sp_objective(model, sum(-item_reward(node)[i] * y[i] for i in 1:5))


# ### Defining the JuDGE model
# The JuDGEModel is defined based on a tree, probability distribution,
# sub_problems, and an optimizer
judy = JuDGEModel(mytree, ConditionallyUniformProbabilities, sub_problems, JuDGE_MP_Solver)

# Solve the model
JuDGE.solve(judy)

# Print the objective, and optimal expansions
println("Objective: "*string(objective_value(judy.master_problem)))
JuDGE.print_expansions(judy,onlynonzero=false)

# Re-solve the subproblems and print the objective
println("Re-solved Objective: " * string(resolve_subproblems(judy)))

# Set up and solve the deterministic equivalent
deteq = DetEqModel(mytree, ConditionallyUniformProbabilities, sub_problems, JuDGE_DE_Solver)
JuDGE.solve(deteq)
println("Deterministic Equivalent Objective: " * string(objective_value(deteq.problem)))
