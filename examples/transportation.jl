using Random, JuMP, JuDGE, Test, DelimitedFiles

include("solvers/setup_gurobi.jl")
#include("solvers/setup_cplex.jl")
#include("solvers/setup_coin.jl")
#include("solvers/setup_glpk.jl")

function transportation()
   mytree = narytree(2,2)
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
      model = Model(JuDGE_SP_Solver)

      @expansion(model, new_supply[supply_nodes]) #invest in more supply
      @expansion(model, new_capacity[supply_nodes,demand_nodes]) #invest in more arc capacity

      @expansioncosts(model, sum(invest_supply_cost(node)[i] * new_supply[i] for i in supply_nodes) +
         sum(invest_arc_cost(node)[i,j] * new_capacity[i,j] for i in supply_nodes for j in demand_nodes))

      @variable(model, x[supply_nodes, demand_nodes] >= 0)

      @sp_objective(model, sum(c_dict[i,j]*x[i,j] for i in supply_nodes, j in demand_nodes))
      @constraint(model, SupplyIncrease[i in supply_nodes],
                  sum(x[i,j] for j in demand_nodes) <= supply(node)[i] + s_dict[i]*new_supply[i])

      @constraint(model, CapacityIncrease[i in supply_nodes, j in demand_nodes],
                  x[i,j] <= c_dict[i,j] + c_dict[i,j]*new_capacity[i,j])

      @constraint(model,demandCon[j in demand_nodes], sum(x[i,j] for i in supply_nodes) == demand(node)[j] )

      return model
   end

   function format_output(s::Symbol,values)
      if s==:new_capacity
         output=Dict{Tuple,Float64}()
         for i in supply_nodes
            for j in demand_nodes
               output[i,j]=values[i,j]*c_dict[i,j]
            end
         end
         return output
      elseif s==:new_supply
         output=Dict{String,Float64}()
         for i in supply_nodes
            output[i]=values[i]*s_dict[i]
         end
         return output
      end
      return nothing
   end

   judy = JuDGEModel(mytree, ConditionallyUniformProbabilities, sub_problems, JuDGE_MP_Solver)
   JuDGE.solve(judy,inttol=10^-9)

   println("\nObjective: "*string(objective_value(judy.master_problem))*"\n")
   JuDGE.print_expansions(judy,format=format_output)

   JuDGE.fix_expansions(judy)
   println("\nRe-solved Objective: " * string(resolve_subproblems(judy)))

   JuDGE.write_solution_to_file(judy,joinpath(@__DIR__,"solution.csv"))

   deteq = DetEqModel(mytree, ConditionallyUniformProbabilities, sub_problems, JuDGE_DE_Solver)
   JuDGE.solve(deteq)
   println("Deterministic Equivalent Objective: " * string(objective_value(deteq.problem)))
   return objective_value(judy.master_problem)
end

@test transportation() ≈ 1081.51 atol = 1e-2
