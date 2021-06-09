using Random
using JuMP
using JuDGE

if !isdefined(@__MODULE__, :JuDGE_MP_Solver)
	# Replace this with another file in `/solvers` as appropriate.
	include("solvers/setup_gurobi.jl")
end


function inventory(;depth=4,degree=4,price_array=nothing,visualise=false,risk=RiskNeutral())
	mytree = narytree(depth,degree)

	if price_array==nothing
		Random.seed!(1000)
		price_array=rand(length(collect(mytree)))
	end

	price=Dict(zip(collect(mytree),price_array))

	function sub_problems(n)
	    sp = JuMP.Model(JuDGE_SP_Solver)

	    @state(sp, -50<=Δstock<=50, state_name=stock, lb=0, ub=200, initial=0)

	    @ongoingcosts(sp, 0.01*stock)

	    @variable(sp, buy>=0)
		@variable(sp, sell>=0)
		#@variable(sp, change, Int)
		#@constraint(sp, Δstock==10*change)
	    @constraint(sp, trade, Δstock == buy - sell)

	    @objective(sp, Min, (buy*1.01)*price[n]-(sell/1.01)*price[n])

	    sp
	end

	model = JuDGEModel(mytree, ConditionallyUniformProbabilities, sub_problems, JuDGE_MP_Solver,check=true, risk=risk)

	JuDGE.solve(model,verbose=1)
	JuDGE.print_expansions(model,onlynonzero=true,inttol=10^-5)
	println("\nRe-solved Objective: " * string(resolve_subproblems(model)))

	if visualise
		solution=JuDGE.solution_to_dictionary(model)
		solution[:prices]=price
		JuDGE.visualize_tree(mytree,solution)
	end
	JuDGE.get_objval(model)
end

if !isdefined(@__MODULE__, :running_tests) || !running_tests
	inventory(visualise=true,risk=RiskNeutral())
	inventory(visualise=true,risk=JuDGE.Risk(0.1,bound=0.0))
end
