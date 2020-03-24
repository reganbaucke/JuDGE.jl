using Random, JuMP, JuDGE, Gurobi, Test

# this test is a test based off of the presentation by andy ages ago
function knapsack_fixed()
   mytree = narytree(2,() -> [Leaf(), Leaf()])
   function invest_cost(node)
      if node == mytree
         180.0
      elseif node == mytree.children[1]
         50.0
      elseif node == mytree.children[2]
         60.0
      elseif node == mytree.children[1].children[1]
         40.0
      elseif node == mytree.children[1].children[2]
         60.0
      elseif node == mytree.children[2].children[1]
         10.0
      elseif node == mytree.children[2].children[2]
         10.0
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
      model = Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0))
      set_silent(model)
      @expansion(model, bag)
      @expansioncosts(model, bag*invest_cost(node))
      @variable(model, y[1:5], Bin)
      @expansionconstraint(model, BagExtension, sum(y[i]*item_volume(node)[i] for i in 1:5) <= 3 + 4 * bag)
      @objective(model, Min, sum(-item_reward(node)[i] * y[i] for i in 1:5))
      return model
   end

   judy = JuDGEModel(mytree, ConditionallyUniformProbabilities, sub_problems, optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0))
   JuDGE.solve(judy)

   println("Objective: "*string(objective_value(judy.master_problem)))
   JuDGE.print_expansions(judy,onlynonzero=false)

   JuDGE.fix_expansions(judy)
   println("Re-solved Objective: " * string(JuDGE.resolve_fixed(judy)))

   return objective_value(judy.master_problem)
end

# this test is a test based off of the test in the old JuDGE
function knapsack_random()
   Random.seed!(1)
   # how many investments?
   numinvest = 2;

   # number of items to pick from in the knapsack?
   numitems = 20

   # size of tree?
   degree = 3
   # height = 5
   height = 3

   totalnodes = Int64((degree^(height + 1) + 1)/(degree-1)) - 1

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

   mytree = narytree(height,() -> [Leaf(),Leaf(),Leaf()])

   nodes = collect(mytree)
   function data(node, input)
      input[findall(x -> x == node, nodes)[1], :]
   end

   function sub_problems(node)
      model = Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0))
      @expansion(model, bag[1:numinvest])
      @expansioncosts(model, sum(data(node,investcost)[i] * bag[i] for i in  1:numinvest))

      @variable(model, y[1:numitems], Bin)
      @expansionconstraint(model, BagExtension ,sum( y[i]*data(node,itemvolume)[i] for i in 1:numitems) <= initialcap + sum(bag[i]*investvol[i] for i in 1:numinvest))
      @objective(model, Min, sum(-data(node,itemcost)[i] * y[i] for i in 1:numitems))
      return model
   end

   judy = JuDGEModel(mytree, ConditionallyUniformProbabilities, sub_problems, optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0))
   JuDGE.solve(judy)

   println("Objective: "*string(objective_value(judy.master_problem)))
   JuDGE.print_expansions(judy)

   JuDGE.fix_expansions(judy)
   println("Re-solved Objective: " * string(JuDGE.resolve_fixed(judy)))

   return objective_value(judy.master_problem)
end

@test knapsack_fixed() ≈ -131.25 atol = 1e-3
@test knapsack_random() ≈ -34.749 atol = 1e-3
