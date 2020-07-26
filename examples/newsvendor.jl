using JuMP, JuDGE, Gurobi, Test

env = Gurobi.Env()

function newsvendor(;cost=5.0,price=8.0,demands=[50.0,60.0,70.0],CVaR=JuDGE.RiskNeutral)
   mytree = narytree(1,length(demands))

   n=Int(round(log(maximum(demands))/log(2)+0.501)) # bits for binary expansion

   demand=Dict{AbstractTree,Float64}()
   demand[get_node(mytree,[1])]=0
   for i in 1:length(demands)
      demand[get_node(mytree,[1,i])]=demands[i]
   end

   function sub_problems(node)
      model = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0))
      @expansion(model, papers_ordered[1:n])
      @expansion(model, papers_arrived[1:n])
      @expansioncosts(model, sum(papers_ordered[i]*cost*2^(i-1) for i in 1:n))
      @variable(model, sales>=0)
      @expansionconstraint(model, maxsale1, sales <= sum(papers_arrived[i]*2^(i-1) for i in 1:n))
      @constraint(model, maxsale2, sales <= demand[node])
      @sp_objective(model, -price*sales)
      return model
   end

   function intertemporal(model,tree)
      for i in 1:n
         for node in collect(tree)
            if node!=tree
               @constraint(model,papers_arrived[node][i]<=papers_ordered[tree][i])
            else
               @constraint(model,papers_arrived[node][i]==0.0)
            end
         end
      end
   end

   judy = JuDGEModel(mytree, ConditionallyUniformProbabilities, sub_problems, optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0),intertemporal=intertemporal,CVaR=CVaR)
   JuDGE.solve(judy)

   println("Objective: "*string(objective_value(judy.master_problem)))

   JuDGE.fix_expansions(judy)
   println("Re-solved Objective: " * string(JuDGE.resolve_fixed(judy)))

   function format_output(s::Symbol,values)
      if s==:papers_ordered
         return sum(values[i]*2^(i-1) for i in 1:n)
      elseif s==:papers_arrived
         return 0.0
      end
      return nothing
   end

   JuDGE.print_expansions(judy,format=format_output)

   deteq = DetEqModel(mytree, ConditionallyUniformProbabilities, sub_problems, optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0),intertemporal=intertemporal,CVaR=CVaR)
   JuDGE.solve(deteq)
   println("Deterministic Equivalent Objective: " * string(objective_value(deteq.problem)))

   JuDGE.print_expansions(deteq,format=format_output)

   objective_value(judy.master_problem)
end

@test newsvendor(cost=5.0,price=8.0,demands=(10,80,100),CVaR=JuDGE.RiskNeutral) ≈ -53.333 atol = 1e-3
@test newsvendor(cost=5.0,price=16.0,demands=(10,80,100),CVaR=JuDGE.RiskNeutral) ≈ -513.333 atol = 1e-3
@test newsvendor(cost=5.0,price=8.0,demands=(10,80,100),CVaR=(0.5,0.5)) ≈ -30.0 atol = 1e-3
@test newsvendor(cost=5.0,price=16.0,demands=(10,80,100),CVaR=(0.5,0.5)) ≈ -320.0 atol = 1e-3
