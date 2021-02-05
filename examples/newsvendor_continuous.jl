using JuMP, JuDGE, Test

include("solvers/setup_gurobi.jl")
#include("solvers/setup_cplex.jl")
#include("solvers/setup_coin.jl")
#include("solvers/setup_glpk.jl")

function newsvendor(;depth=1,cost=5.0,price=8.0,demands=[50.0,60.0,70.0],CVaR=JuDGE.RiskNeutral)
   mytree = narytree(depth,length(demands))

   demand=Dict{AbstractTree,Float64}()
   nodes=collect(mytree,order=:breadth)
   for i in eachindex(nodes)
      if i==1
         demand[nodes[i]]=0
      else
         demand[nodes[i]]=demands[(i-2)%length(demands)+1]
      end
   end

   function sub_problems(node)
      model = Model(JuDGE_SP_Solver)
      JuDGE.@expansion_cont(model, 0<=papers_ordered<=1000, 1, 1)
      @capitalcosts(model, papers_ordered*cost)
      @variable(model, sales>=0)
      @constraint(model, maxsale1, sales <= papers_ordered)
      @constraint(model, maxsale2, sales <= demand[node])
      @objective(model, Min, -price*sales)
      return model
   end

   judy = JuDGEModel(mytree, ConditionallyUniformProbabilities, sub_problems, JuDGE_MP_Solver, CVaR=CVaR)
   JuDGE.solve(judy)

   println("Objective: "*string(objective_value(judy.master_problem)))

   println("Re-solved Objective: " * string(resolve_subproblems(judy)))

   JuDGE.print_expansions(judy)
   JuDGE.write_solution_to_file(judy,joinpath(@__DIR__,"solution.csv"))
   deteq = DetEqModel(mytree, ConditionallyUniformProbabilities, sub_problems, JuDGE_DE_Solver, CVaR=CVaR)
   JuDGE.solve(deteq)
   println("Deterministic Equivalent Objective: " * string(objective_value(deteq.problem)))

   JuDGE.print_expansions(deteq)
   JuDGE.write_solution_to_file(deteq,joinpath(@__DIR__,"deteq_solution.csv"))
   objective_value(judy.master_problem)
end

@test newsvendor(cost=5.0,price=8.0,demands=[10,80,100],CVaR=JuDGE.RiskNeutral) ≈ -53.333 atol = 1e-3
@test newsvendor(cost=5.0,price=16.0,demands=[10,80,100],CVaR=JuDGE.RiskNeutral) ≈ -513.333 atol = 1e-3
@test newsvendor(cost=5.0,price=8.0,demands=[10,80,100],CVaR=(0.5,0.5)) ≈ -30.0 atol = 1e-3
@test newsvendor(cost=5.0,price=16.0,demands=[10,80,100],CVaR=(0.5,0.5)) ≈ -320.0 atol = 1e-3
@test newsvendor(depth=2,cost=5.0,price=8.0,demands=[10,20,30],CVaR=(0.15,0.5)) ≈ -61.028 atol = 1e-3
@test newsvendor(depth=3,cost=5.0,price=8.0,demands=[10,20,30],CVaR=(0.15,0.05)) ≈ -90.617 atol = 1e-3
