using Random, JuMP, JuDGE, Gurobi, Test

env = Gurobi.Env()
# this test is a test based off of the presentation by andy ages ago
function knapsack_fixed()
   mytree = narytree(2,2)

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

   ### with judge
   function sub_problems(node)
      model = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0))
      @expansion(model, bag)
      @expansioncosts(model, bag*invest_cost(node))
      @variable(model, y[1:5], Bin)
      @expansionconstraint(model, BagExtension, sum(y[i]*item_volume(node)[i] for i in 1:5) <= 3 + 4 * bag)
      @sp_objective(model, sum(-item_reward(node)[i] * y[i] for i in 1:5))
      return model
   end

   judy = JuDGEModel(mytree, ConditionallyUniformProbabilities, sub_problems, optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0))
   JuDGE.solve(judy)

   println("Objective: "*string(objective_value(judy.master_problem)))
   JuDGE.print_expansions(judy,onlynonzero=false)

   JuDGE.fix_expansions(judy)
   println("Re-solved Objective: " * string(JuDGE.resolve_fixed(judy)))

   deteq = DetEqModel(mytree, ConditionallyUniformProbabilities, sub_problems, optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0))
   JuDGE.solve(deteq)
   println("Deterministic Equivalent Objective: " * string(objective_value(deteq.problem)))

   return objective_value(judy.master_problem)
end

# this test is a test based off of the test in the old JuDGE
# this test is a test based off of the test in the old JuDGE
function knapsack_random()
   Random.seed!(1)
   # how many investments?
   numinvest = 2;

   # number of items to pick from in the knapsack?
   numitems = 20

   # size of tree?
   degree = 3
   height = 3

   totalnodes = Int64((degree^(height+1) - 1)/(degree-1))

   investcost = zeros(totalnodes,numinvest)
   for i = 1:totalnodes
      investcost[i,:] = (rand(numinvest)*2  + 2*[2.0,3.5])*(1-((i-1)/(totalnodes*1.2)))
   end

   # investvol = [40,45,50,70]
   investvol = [40,50]
   initialcap = 80

   itemvolume = zeros(totalnodes,numitems)
   for i = 1:totalnodes
      itemvolume[i,:] = ((rand(numitems) .- 0.5)*2)*2 + collect(range(4,22,length = numitems))
   end

   itemcost = zeros(totalnodes,numitems)
   for i = 1:totalnodes
      itemcost[i,:] = ((rand(numitems) .- 0.5)*2)*0.5 + collect(range(0.5,1,length = numitems))
   end

   mytree = narytree(height,degree)

   nodes = collect(mytree)
   function data(node, input)
      input[findall(x -> x == node, nodes)[1], :]
   end

   function sub_problems(node)
      model = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0))
      @expansion(model, bag[1:numinvest])
      @expansioncosts(model, sum(data(node,investcost)[i] * bag[i] for i in  1:numinvest))
      @variable(model, y[1:numitems], Bin)
      @expansionconstraint(model, BagExtension ,sum( y[i]*data(node,itemvolume)[i] for i in 1:numitems) <= initialcap + sum(bag[i]*investvol[i] for i in 1:numinvest))
      @sp_objective(model, sum(-data(node,itemcost)[i] * y[i] for i in 1:numitems))
      return model
   end

   judy = JuDGEModel(mytree, ConditionallyUniformProbabilities, sub_problems, optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0))
   JuDGE.solve(judy)

   println("Objective: "*string(objective_value(judy.master_problem)))
   JuDGE.print_expansions(judy)

   JuDGE.fix_expansions(judy)
   println("Re-solved Objective: " * string(JuDGE.resolve_fixed(judy)))

   deteq = DetEqModel(mytree, ConditionallyUniformProbabilities, sub_problems, optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0))
   JuDGE.solve(deteq)
   println("Deterministic Equivalent Objective: " * string(objective_value(deteq.problem)))

   return objective_value(judy.master_problem)
end

function knapsack_branch_and_price()
   Random.seed!(50)
   # how many investments?
   numinvest = 5;

   # number of items to pick from in the knapsack?
   numitems = 10

   # size of tree?
   degree = 3
   height = 4

   totalnodes = Int64((degree^(height+1) - 1)/(degree-1))

   investcost = zeros(totalnodes,numinvest)
   for i = 1:totalnodes
      investcost[i,:] = ([1,1.8,3.5,6.8,13.5])*(1-((i-1)/(totalnodes*1.2)))
   end

   # investvol = [40,45,50,70]
   investvol = [1,2,4,8,16]
   initialcap = 0

   itemvolume = zeros(totalnodes,numitems)
   for i = 1:totalnodes
      itemvolume[i,:] = ((rand(numitems))*2) + collect(range(4,22,length = numitems))
   end

   itemcost = zeros(totalnodes,numitems)
   for i = 1:totalnodes
      itemcost[i,:] = ((rand(numitems) .- 0.5)*2)*2# + collect(range(0.5,1,length = numitems))
   end

   mytree = narytree(height,degree)

   nodes = collect(mytree)
   function data(node, input)
      input[findall(x -> x == node, nodes)[1], :]
   end

   function sub_problems(node)
      model = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0,"MIPGap" => 0.0))
      @expansion(model, bag[1:numinvest])
      @expansioncosts(model, sum(data(node,investcost)[i] * bag[i] for i in  1:numinvest))
      @variable(model, y[1:numitems], Bin)
      @expansionconstraint(model, BagExtension ,sum( y[i]*data(node,itemvolume)[i] for i in 1:numitems) <= initialcap + sum(bag[i]*investvol[i] for i in 1:numinvest))
      @sp_objective(model, sum(-data(node,itemcost)[i] * y[i] for i in 1:numitems))
      return model
   end

   judy = JuDGEModel(mytree, ConditionallyUniformProbabilities, sub_problems,
                     optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0, "Method" => 2, "Crossover" => 0,  "MIPGap" => 0.0))

   best=JuDGE.branch_and_price(judy,rlx_abstol=10^-7,inttol=10^-6,
   					branch_method=JuDGE.constraint_branch,search=:lowestLB)

   println("Objective: "*string(objective_value(best.master_problem)))
   JuDGE.print_expansions(best)

   deteq = DetEqModel(mytree, ConditionallyUniformProbabilities, sub_problems, optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0))
   JuDGE.solve(deteq)
   println("Deterministic Equivalent Objective: " * string(objective_value(deteq.problem)))

   return objective_value(best.master_problem)
end

function knapsack_risk_averse()
   Random.seed!(50)
   # how many investments?
   numinvest = 5;

   # number of items to pick from in the knapsack?
   numitems = 10

   # size of tree?
   degree = 3
   height = 4

   totalnodes = Int64((degree^(height+1) - 1)/(degree-1))

   investcost = zeros(totalnodes,numinvest)
   for i = 1:totalnodes
      investcost[i,:] = ([1,1.8,3.5,6.8,13.5])*(1-((i-1)/(totalnodes*1.2)))
   end

   # investvol = [40,45,50,70]
   investvol = [1,2,4,8,16]
   initialcap = 0

   itemvolume = zeros(totalnodes,numitems)
   for i = 1:totalnodes
      itemvolume[i,:] = ((rand(numitems))*2) + collect(range(4,22,length = numitems))
   end

   itemcost = zeros(totalnodes,numitems)
   for i = 1:totalnodes
      itemcost[i,:] = ((rand(numitems) .- 0.5)*2)*2# + collect(range(0.5,1,length = numitems))
   end

   mytree = narytree(height,degree)

   nodes = collect(mytree)
   function data(node, input)
      input[findall(x -> x == node, nodes)[1], :]
   end

   function sub_problems(node)
      model = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0,"MIPGap" => 0.0))
      @expansion(model, bag[1:numinvest])
      @expansioncosts(model, sum(data(node,investcost)[i] * bag[i] for i in  1:numinvest))
      @variable(model, y[1:numitems], Bin)
      @expansionconstraint(model, BagExtension ,sum( y[i]*data(node,itemvolume)[i] for i in 1:numitems) <= initialcap + sum(bag[i]*investvol[i] for i in 1:numinvest))
      @sp_objective(model, sum(-data(node,itemcost)[i] * y[i] for i in 1:numitems))
      return model
   end

   judy = JuDGEModel(mytree, ConditionallyUniformProbabilities, sub_problems,
                     optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0, "Method" => 2, "Crossover" => 0,  "MIPGap" => 0.0),CVaR=(0.5,0.05))

   best=JuDGE.branch_and_price(judy,rlx_abstol=10^-6,inttol=10^-6,
   					branch_method=JuDGE.constraint_branch,search=:lowestLB)

   println("Objective: "*string(objective_value(best.master_problem)))
   JuDGE.print_expansions(best)

   JuDGE.fix_expansions(best)
   println("Re-solved Objective: " * string(JuDGE.resolve_fixed(best)))
   deteq = DetEqModel(mytree, ConditionallyUniformProbabilities, sub_problems, optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0, "MIPGap" => 0.0),CVaR=(0.5,0.05))
   JuDGE.solve(deteq)
   println("Deterministic Equivalent Objective: " * string(objective_value(deteq.problem)))

   return objective_value(best.master_problem)
end

function knapsack_delayed_investment(;CVaR=(0.0,1.0))
   Random.seed!(100)
   # how many investments?
   numinvest = 2;

   # number of items to pick from in the knapsack?
   numitems = 25

   # size of tree?
   degree = 3
   height = 3

   totalnodes = Int64((degree^(height+1) - 1)/(degree-1))

   investcost = zeros(totalnodes,numinvest)
   for i = 1:totalnodes
      investcost[i,:] = (rand(numinvest)*2  + 2*[2.0,3.5])*(1-((i)/(totalnodes)))
   end

   # investvol = [40,45,50,70]
   investvol = [40,50]
   initialcap = 80

   itemvolume = zeros(totalnodes,numitems)
   for i = 1:totalnodes
      itemvolume[i,:] = ((rand(numitems) .- 0.5)*2)*2 + collect(range(4,22,length = numitems))
   end

   itemcost = zeros(totalnodes,numitems)
   for i = 1:totalnodes
      itemcost[i,:] = ((rand(numitems) .- 0.5)*2)*0.5 + collect(range(0.5,1,length = numitems))
   end

   mytree = narytree(height,degree)

   nodes = collect(mytree)
   function data(node, input)
      input[findall(x -> x == node, nodes)[1], :]
   end

   function sub_problems(node)
      model = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0))
      @expansion(model, bag_bought[1:numinvest])
      @expansion(model, bag_arrived[1:numinvest])
      @expansioncosts(model, sum(data(node,investcost)[i] * bag_bought[i] for i in  1:numinvest))
      @variable(model, y[1:numitems], Bin)
      @expansionconstraint(model, BagExtension ,sum( y[i]*data(node,itemvolume)[i] for i in 1:numitems) <= initialcap + sum(bag_arrived[i]*investvol[i] for i in 1:numinvest))
      @sp_objective(model, sum(-data(node,itemcost)[i] * y[i] for i in 1:numitems))
      return model
   end

   function intertemporal(model,tree,node,current_expansions,previous_expansions)
      if previous_expansions==nothing
         for i in eachindex(current_expansions[:bag_arrived])
            @constraint(model,current_expansions[:bag_arrived][i]==0)
         end
      else
         #parent_fn=JuDGE.parent_builder(tree)
         #p=parent_fn(node)
         for i in eachindex(current_expansions[:bag_arrived])
            @constraint(model,current_expansions[:bag_arrived][i]<=sum(prev[:bag_bought][i] for (n,prev) in previous_expansions))
         end
      end
   end

   judy = JuDGEModel(mytree, ConditionallyUniformProbabilities, sub_problems, optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0),intertemporal=intertemporal,CVaR=CVaR)
   #JuDGE.solve(judy)
   judy=JuDGE.branch_and_price(judy,rlx_abstol=10^-6,inttol=10^-6,
                  branch_method=JuDGE.constraint_branch,search=:lowestLB)
   println("Objective: "*string(objective_value(judy.master_problem)))
   JuDGE.print_expansions(judy)

   JuDGE.fix_expansions(judy)
   println("Re-solved Objective: " * string(JuDGE.resolve_fixed(judy)))

   deteq = DetEqModel(mytree, ConditionallyUniformProbabilities, sub_problems, optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0, "MIPGap" => 0.0),intertemporal=intertemporal,CVaR=CVaR)
   JuDGE.solve(deteq)
   println("Deterministic Equivalent Objective: " * string(objective_value(deteq.problem)))
   return objective_value(judy.master_problem)
end

function knapsack_divestment()
   mytree = narytree(2,2)
   function divest_revenue(node)
      if node == mytree
         60.0
      elseif node == mytree.children[1]
         58.0
      elseif node == mytree.children[2]
         35.0
      elseif node == mytree.children[1].children[1]
         40.0
      elseif node == mytree.children[1].children[2]
         35.0
      elseif node == mytree.children[2].children[1]
         25.0
      elseif node == mytree.children[2].children[2]
         15.0
      end
   end

   function item_volume(node)
      if node == mytree
         [6, 2, 1, 1, 1]
      elseif node == mytree.children[1]
         [8, 2, 2, 2, 1]
      elseif node == mytree.children[2]
         [8, 1, 1, 1, 3]
      elseif node == mytree.children[1].children[1]
         [4, 4, 3, 1, 2]
      elseif node == mytree.children[1].children[2]
         [1, 3, 1, 1, 2]
      elseif node == mytree.children[2].children[1]
         [7, 3, 1, 1, 1]
      elseif node == mytree.children[2].children[2]
         [2, 5, 2, 1, 2]
      end
   end

   function item_reward(node)
      if node == mytree
         [60, 20, 10, 15, 10]
      elseif node == mytree.children[1]
         [8, 10, 20, 20, 10]
      elseif node == mytree.children[2]
         [8, 10, 15, 10, 30]
      elseif node == mytree.children[1].children[1]
         [40, 40, 35, 10, 20]
      elseif node == mytree.children[1].children[2]
         [15, 35, 15, 15, 20]
      elseif node == mytree.children[2].children[1]
         [70, 30, 15, 15, 10]
      elseif node == mytree.children[2].children[2]
         [25, 50, 25, 15, 20]
      end
   end

   ### with judge
   function sub_problems(node)
      model = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0))
      set_silent(model)
      @forced_expansion(model, bag)
      @expansioncosts(model, -bag*divest_revenue(node))
      @variable(model, y[1:5], Bin)
      @expansionconstraint(model, BagExtension, sum(y[i]*item_volume(node)[i] for i in 1:5) <= 4 - 2 * bag)
      @sp_objective(model, sum(-item_reward(node)[i] * y[i] for i in 1:5))
      return model
   end

   judy = JuDGEModel(mytree, ConditionallyUniformProbabilities, sub_problems, optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0))
   JuDGE.solve(judy)

   println("Objective: "*string(objective_value(judy.master_problem)))
   JuDGE.print_expansions(judy,onlynonzero=false)

   JuDGE.fix_expansions(judy)
   println("Re-solved Objective: " * string(JuDGE.resolve_fixed(judy)))

   deteq = DetEqModel(mytree, ConditionallyUniformProbabilities, sub_problems, optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0))
   JuDGE.solve(deteq)
   println("Deterministic Equivalent Objective: " * string(objective_value(deteq.problem)))

   return objective_value(judy.master_problem)
end

@test knapsack_fixed() ≈ -131.25 atol = 1e-3
@test knapsack_random() ≈ -34.749 atol = 1e-3
@test knapsack_branch_and_price() ≈ -0.69456 atol = 1e-4
@test knapsack_risk_averse() ≈ -0.27292 atol = 1e-4
@test knapsack_delayed_investment(CVaR=(0.0,1.0)) ≈ -34.058 atol = 1e-3
@test knapsack_delayed_investment(CVaR=(0.95,0.05)) ≈ -31.344 atol = 1e-3
@test knapsack_divestment() ≈ -145.25 atol = 1e-3
