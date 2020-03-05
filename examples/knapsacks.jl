using Random
using JuMP
using JuDGE
using Gurobi

# one node version of the problem
function test_one()
   items = 1:5
   itemreward = Dict(i => 1.5*i for i in items)
   itemvolume = Dict(i => 10*i for i in items)
   initialcap = 50

   expansionvolume = 15
   expansioncost = 0

   # the deterministic equivalent
   det = Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0))

   @variable(det, x[i in items], Bin)

   @variable(det, extensionbag, Bin)
   @constraint(det, sum( x[i]*itemvolume[i] for i in items) <= initialcap + extensionbag*expansionvolume)
   @objective(det, Min, sum(-itemreward[i]*x[i] for i in items) + expansioncost*extensionbag )

   optimize!(det)

   function blob()
      ## JuDGE Version
      mytree = Leaf()
      probabilities = x -> (y -> 1.0)
      function sub_problem_builder(node)
         sub = Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0))
         @expansion(sub, extensionbag)
         @expansioncosts(sub, expansioncost * extensionbag)

         @variable(sub, x[i in items], Bin)
         @expansionconstraint(sub, BagExtension , sum( x[i]*itemvolume[i] for i in items) <= initialcap + extensionbag*expansionvolume)
         @objective(sub, Min, sum(-itemreward[i]*x[i] for i in items))
         return sub
      end

      hello = JuDGEModel(mytree, probabilities, sub_problem_builder, optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0))

      judgesolve(hello)
      return hello
   end
   hello = blob()

   println("-----------")
   println(objective_value(det))
   println(objective_value(hello.master_problem))
   println("-----------")
end

# this test is a test based off of the test in the old JuDGE
function test_two()

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

   judgesolve(judy)
end

# this test is a test based off of the presentation by andy ages ago
function test_three()
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
    judgesolve(judy)
    objective_value(judy.master_problem)
   judy
   ###

   ### without judge

   # model = Model(with_optimizer(Gurobi.Optimizer))
   # @variable(model, y[1:5, nodes in collect(mytree)], Bin)
   # @variable(model, bag[nodes in collect(mytree)], Bin)

   # @constraint(model, [node in collect(mytree)] ,sum(bag[past] for past in history(mytree)(node)) <= 1)
   # @constraint(model, [node in collect(mytree)], sum(y[i,node]*item_volume(node)[i] for i in 1:5) <= 3 + 4 * sum(bag[past] for past in history(mytree)(node)))

   # @objective(model,Min, sum(ConditionallyUniformProbabilities(mytree)(node)*(sum(-item_reward(node)[i]*y[i,node] for i in 1:5) + bag[node]*invest_cost(node)) for node in collect(mytree)))

   # (model,judy)

end
