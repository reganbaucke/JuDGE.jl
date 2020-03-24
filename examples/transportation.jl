using Random, JuMP, JuDGE, Gurobi, Test, DelimitedFiles

function transportation()
   mytree = narytree(2,() -> [Leaf(), Leaf()])
   get_parent=JuDGE.parent_builder(mytree)

   function invest_supply_cost(node)
      if get_parent(node)==nothing
         return Dict( zip( supply_nodes, [1.0,2.0]) )
      else
         p=get_parent(node)
         for i in 1:length(p.children)
            if p.children[i]==node
               temp=deepcopy(invest_supply_cost(p))
               for key in keys(temp)
                  temp[key]*=(0.2*i+0.8)
               end
               return temp
            end
         end
      end
   end

   function invest_arc_cost(node)
      if get_parent(node)==nothing
         temp=[]
         for i in supply_nodes
            for j in demand_nodes
               push!(temp,(i,j))
            end
         end
         return Dict( zip( temp, [1.0,2.0,3.0,4.0,5.0,6.0]) )
      else
         p=get_parent(node)
         for i in 1:length(p.children)
            if p.children[i]==node
               temp=deepcopy(invest_arc_cost(p))
               for key in keys(temp)
                  temp[key]*=(0.2*i+0.8)
               end
               return temp
            end
         end
      end
   end

   function demand(node)
      if get_parent(node)==nothing
         return d_dict
      else
         p=get_parent(node)
         for i in 1:length(p.children)
            if p.children[i]==node
               temp=deepcopy(demand(p))
               for key in keys(temp)
                  temp[key]+=i*i
               end
               return temp
            end
         end
      end
   end

   function supply(node)
      return s_dict
   end

   data_file = "transportation.csv"
   data = readdlm(joinpath(@__DIR__,data_file),',')

   supply_nodes = data[3:end, 2]
   s = data[3:end, 1]

   demand_nodes = collect(data[2, 3:end])
   d = collect(data[1, 3:end])

   c = data[3:end, 3:end]

   # Converting arrays to dictionaries
   s_dict = Dict( zip( supply_nodes, s*2) )
   d_dict = Dict( zip( demand_nodes, d) )

   c_dict = Dict()
   for i in 1:length(supply_nodes)
     for j in 1:length(demand_nodes)
       c_dict[supply_nodes[i], demand_nodes[j]] = c[i,j]
     end
   end

   ### with judge
   function sub_problems(node)
      model = Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0))

      @expansion(model, new_supply[supply_nodes]) #invest in more supply
      @expansion(model, new_capacity[supply_nodes,demand_nodes]) #invest in more arc capacity

      @expansioncosts(model, sum(invest_supply_cost(node)[i] * new_supply[i] for i in supply_nodes) +
         sum(invest_arc_cost(node)[i,j] * new_capacity[i,j] for i in supply_nodes for j in demand_nodes))

      @variable(model, x[supply_nodes, demand_nodes] >= 0)

      @objective(model, Min, sum(c_dict[i,j]*x[i,j] for i in supply_nodes, j in demand_nodes))
      @expansionconstraint(model, SupplyIncrease[i in supply_nodes],
                  sum(x[i,j] for j in demand_nodes) <= supply(node)[i] + s_dict[i]*new_supply[i])

      @expansionconstraint(model, CapacityIncrease[i in supply_nodes, j in demand_nodes],
                  x[i,j] <= c_dict[i,j] + c_dict[i,j]*new_capacity[i,j])

      @constraint(model,demandCon[j in demand_nodes], sum(x[i,j] for i in supply_nodes) == demand(node)[j] )

      return model
   end

   judy = JuDGEModel(mytree, ConditionallyUniformProbabilities, sub_problems, optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" => 0))
   JuDGE.solve(judy)

   println("\nObjective: "*string(objective_value(judy.master_problem))*"\n")
   JuDGE.print_expansions(judy)

   JuDGE.fix_expansions(judy)
   println("\nRe-solved Objective: " * string(JuDGE.resolve_fixed(judy)))

   JuDGE.write_solution_to_file(judy,joinpath(@__DIR__,"solution.csv"))

   return objective_value(judy.master_problem)
end

@test transportation() ≈ 1081.51 atol = 1e-2
