using JuMP, JuDGE, Test

include("solvers/setup_gurobi.jl")
#include("solvers/setup_cplex.jl")
#include("solvers/setup_coin.jl")
#include("solvers/setup_glpk.jl")

function newsvendor(;depth=1,cost=5.0,price=8.0,demands=[50.0,60.0,70.0],CVaR=JuDGE.RiskNeutral,verbose=2)
   mytree = narytree(depth,length(demands))

   n=Int(round(log(maximum(demands))/log(2)+0.501)) # bits for binary expansion

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
      @expansion(model, papers_ordered[1:n], Bin, lag=1, duration=1)
      @capitalcosts(model, sum(papers_ordered[i]*cost*2^(i-1) for i in 1:n))
      @variable(model, sales>=0)
      @constraint(model, maxsale1, sales <= sum(papers_ordered[i]*2^(i-1) for i in 1:n))
      @constraint(model, maxsale2, sales <= demand[node])
      @objective(model, Min, -price*sales)
      return model
   end

   judy = JuDGEModel(mytree, ConditionallyUniformProbabilities, sub_problems, JuDGE_MP_Solver, risk=CVaR)
   judy=JuDGE.branch_and_price(judy,search=:lowestLB,branch_method=JuDGE.variable_branch,verbose=verbose)

   println("Objective: "*string(objective_value(judy.master_problem)))

   println("Re-solved Objective: " * string(resolve_subproblems(judy)))

   function format_output(s::Symbol,values)
      if s==:papers_ordered
         return sum(values[i]*2^(i-1) for i in 1:n)
      end
      return nothing
   end

   JuDGE.print_expansions(judy,format=format_output)
   JuDGE.write_solution_to_file(judy,joinpath(@__DIR__,"solution.csv"))
   deteq = DetEqModel(mytree, ConditionallyUniformProbabilities, sub_problems, JuDGE_DE_Solver, risk=CVaR)
   JuDGE.solve(deteq)
   println("Deterministic Equivalent Objective: " * string(objective_value(deteq.problem)))

   JuDGE.print_expansions(deteq,format=format_output)
   JuDGE.write_solution_to_file(deteq,joinpath(@__DIR__,"deteq_solution.csv"))
   objective_value(judy.master_problem)
end

@test newsvendor(cost=5.0,price=8.0,demands=[10,80,100],CVaR=JuDGE.RiskNeutral) ≈ -53.333 atol = 1e-3
@test newsvendor(cost=5.0,price=16.0,demands=[10,80,100],CVaR=JuDGE.RiskNeutral) ≈ -513.333 atol = 1e-3
@test newsvendor(cost=5.0,price=8.0,demands=[10,80,100],CVaR=(0.5,0.5)) ≈ -30.0 atol = 1e-3
@test newsvendor(cost=5.0,price=16.0,demands=[10,80,100],CVaR=(0.5,0.5)) ≈ -320.0 atol = 1e-3
@test newsvendor(depth=2,cost=5.0,price=8.0,demands=[10,20,30],CVaR=(0.15,0.5),verbose=1) ≈ -61.011 atol = 1e-3
@test newsvendor(depth=3,cost=5.0,price=8.0,demands=[10,20,30],CVaR=(0.15,0.05),verbose=0) ≈ -90.526 atol = 1e-3
