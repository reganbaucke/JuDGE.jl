module JuDGE

using MathOptInterface
using Printf
using Distributed
@everywhere using JuMP

include("tree.jl")
include("macros.jl")
include("convergence.jl")
include("utilities.jl")
include("deteq.jl")

mutable struct Bounds
	UB::Float64
	LB::Float64
end

struct Column
	node::AbstractTree
	singleVars::Array{Any,1}
	arrayVars::Dict{Any,Array{Any,1}}
	obj::Float64
end

const RiskNeutral=(0.0,1.0)

struct JuDGEModel
	tree::AbstractTree
	master_problem::JuMP.Model
	sub_problems::Dict{AbstractTree,JuMP.Model}
	bounds::Bounds
	discount_factor::Float64
	master_solver
	probabilities
	CVaR::Tuple{Float64,Float64}
	sideconstraints
	log::Array{ConvergenceState,1}
	optimizer_settings::Array{Int64,1}
end

function JuDGEModel(tree::T where T <: AbstractTree, probabilities::Dict{AbstractTree,Float64}, sub_problems::Dict{AbstractTree,JuMP.Model}, solver, bounds::Bounds, discount_factor::Float64, CVaR::Tuple{Float64,Float64}, sideconstraints, optimizer_settings::Array{Int64,1})
	master_problem = build_master(sub_problems, tree, probabilities, solver, discount_factor, CVaR, sideconstraints)
	JuDGEModel(tree,master_problem,sub_problems,bounds,discount_factor,solver,probabilities,CVaR,sideconstraints,Array{ConvergenceState,1}(),optimizer_settings)
end

"""
	JuDGEModel(tree::AbstractTree,
               probabilities,
               sub_problem_builder::Function,
               solver;
               discount_factor=1.0,
               CVaR=RiskNeutral,
               sideconstraints=nothing,
			   parallel=false,
			   sp_solver=nothing,
			   check=true
			   )

Define a JuDGE model.

### Required arguments
`tree` is a reference to a scenario tree

`probabilities` is either a function, which returns a dictionary of the probabilities
of all nodes in a tree, or simply the dictionary itself

`sub_problem_builder` is a function mapping a node to a JuMP model for each subproblems

`solver` is a reference to the optimizer used for the master problem (with appropriate settings);
 this can also be a tuple containing two optimizers (one for solving the relaxation, and one for
 solving the binary model)

### Optional arguments
`discount_factor` is a number between 0 and 1 defining a constant discount factor along each arc
in the scenario tree

`CVaR` is a tuple with the two CVaR parameters: (λ, β)

`sideconstraints` is a function which specifies side constraints in the master problem, see
[Tutorial 9: Side-constraints](@ref) for further details

`parallel` is a boolean, setting whether the sub-problems will be formulated in parallel

`sp_solver` if formulating the sub-problems in parallel, some solvers have issues; in these case
the sub-problem solver can be set using this argument

`check` is a boolean, which can be set to `false` to disable the validation of the JuDGE model.

### Examples
	judge = JuDGEModel(tree, ConditionallyUniformProbabilities, sub_problems,
                                    Gurobi.Optimizer)
	judge = JuDGEModel(tree, probabilities, sub_problems, CPLEX.Optimizer,
                                    discount_factor=0.9, CVaR=(0.5,0.1)))
"""
function JuDGEModel(tree::T where T <: AbstractTree, probabilities, sub_problem_builder::Function, solver; discount_factor=1.0, CVaR=(0.0,1.0), sideconstraints=nothing, parallel=false, sp_solver=nothing, check=true)
	println("")
	println("Establishing JuDGE model for tree: " * string(tree))
	if typeof(probabilities) <: Function
		probabilities = probabilities(tree)
	end
	if typeof(probabilities)!=Dict{AbstractTree,Float64}
		error("\'probabilities\' needs to be a dictionary mapping AbstractTree to Float64\nor a function that generates such a dictionary")
	end

	nodes=collect(tree)
	if parallel
		sps=pmap(sub_problem_builder,nodes)
		i=1
		sub_problems=Dict()
		for n in nodes
			sub_problems[n]=sps[i]
			i+=1
		end
	else
		sub_problems = Dict(i => sub_problem_builder(i) for i in nodes)
	end

	if sp_solver!=nothing
		for n in nodes
			set_optimizer(sub_problems[n],sp_solver)
		end
	end

	if check
		print("Checking sub-problem format...")
		check_specification_is_legal(sub_problems)
		println("Passed")
	else
		println("Skipping checks of sub-problem format")
	end
	scale_objectives(tree,sub_problems,discount_factor)
	print("Building master problem...")
	master_problem = build_master(sub_problems, tree, probabilities, solver, discount_factor, CVaR, sideconstraints)
	println("Complete")
	JuDGEModel(tree,master_problem,sub_problems,Bounds(Inf,-Inf),discount_factor,solver,probabilities,CVaR,sideconstraints,Array{ConvergenceState,1}(),Int64[])
end

include("branchandprice.jl")

function build_master(sub_problems::T where T <: Dict, tree::T where T <: AbstractTree, probabilities::Dict{AbstractTree,Float64}, solver, discount_factor::Float64, CVaR::Tuple{Float64,Float64}, sideconstraints)
	if typeof(solver) <: Tuple
		model = Model(solver[1])
	else
		model = Model(solver)
	end
	@objective(model,Min,0)

	history_function = history(tree)
	depth_function = depth(tree)

	model.ext[:columns] = Array{Column,1}()
	leafs=Array{Leaf,1}()

	model.ext[:expansions] = Dict{AbstractTree,Dict{Symbol,Any}}()

	for (node,sp) in sub_problems
		if typeof(node)==Leaf
			push!(leafs,node)
		end
		model.ext[:expansions][node] = Dict{Symbol,Any}()
		for (name,variable) in sp.ext[:expansions]
			model.ext[:expansions][node][name] = copy_variable!(model, variable, relaxbinary)
		end
	end

	model.ext[:scenprofit_var] = Dict{Leaf,VariableRef}()
	model.ext[:scenprofit_con] = Dict{Leaf,ConstraintRef}()

	for leaf in leafs
		model.ext[:scenprofit_var][leaf] = @variable(model)
		model.ext[:scenprofit_con][leaf] = @constraint(model, 0 == model.ext[:scenprofit_var][leaf])
		set_objective_coefficient(model, model.ext[:scenprofit_var][leaf], probabilities[leaf])
	end

    if CVaR[1]>0.0 && CVaR[2]<1.0
		eta=@variable(model)
	    for leaf in leafs
			v = @variable(model)
			w = @variable(model)
			set_lower_bound(v,0.0)
			set_lower_bound(w,0.0)
			@constraint(model,v>=eta-model.ext[:scenprofit_var][leaf])
			@constraint(model,w>=model.ext[:scenprofit_var][leaf]-eta)
			set_objective_coefficient(model, v, probabilities[leaf]*CVaR[1])
			set_objective_coefficient(model, w, probabilities[leaf]*CVaR[1]/CVaR[2]*(1-CVaR[2]))
		end
	end

	for leaf in leafs
		nodes=history_function(leaf)
		for n in eachindex(nodes)
			node=nodes[n]
			sp=sub_problems[node]
			df=discount_factor^depth_function(node)
			for (name,variable) in sp.ext[:expansions]
				interval=max(1,n-sp.ext[:options][name][3]-sp.ext[:options][name][2]+1):n-sp.ext[:options][name][2]
				disc=Dict{Int64,Float64}()
				for i in interval
					disc[i]=df/discount_factor^(i-1)
				end
				if typeof(variable) <: AbstractArray
					for i in eachindex(variable)
						cost_coef=df*coef(sp.ext[:capitalcosts],variable[i])
						for j in interval
							cost_coef+=disc[j]*coef(sp.ext[:ongoingcosts],variable[i])
						end
						set_normalized_coefficient(model.ext[:scenprofit_con][leaf], model.ext[:expansions][node][name][i], cost_coef)
					end
				else
					cost_coef=df*coef(sp.ext[:capitalcosts],variable)
					for j in interval
						cost_coef+=disc[j]*coef(sp.ext[:ongoingcosts],variable)
					end
					set_normalized_coefficient(model.ext[:scenprofit_con][leaf], model.ext[:expansions][node][name], cost_coef)
				end
			end
		end
	end

	# create the cover constraints
	model.ext[:coverconstraint] = Dict{AbstractTree,Dict{Symbol,Any}}()
	for (node,sp) in sub_problems
		model.ext[:coverconstraint][node] = Dict{Symbol,Any}()
		past=history_function(node)
		for (name,variable) in sp.ext[:expansions]
			interval=sp.ext[:options][name][2]+1:min(sp.ext[:options][name][2]+sp.ext[:options][name][3],length(past))
			if typeof(variable) <: AbstractArray
				model.ext[:coverconstraint][node][name] = Dict()
				for i in eachindex(variable)
					if sp.ext[:options][name][1]
						model.ext[:coverconstraint][node][name][i] = @constraint(model, 0 >= sum(model.ext[:expansions][past[index]][name][i] for index in interval))
					else
						model.ext[:coverconstraint][node][name][i] = @constraint(model, 0 <= sum(model.ext[:expansions][past[index]][name][i] for index in interval))
					end
				end
				# if typeof(node)==Leaf && sp.ext[:options][name][1]
				# 	for i in eachindex(variable)
				# 		@constraint(model, sum(model.ext[:expansions][past][name][i] for past in history_function(node))<=1)
				# 	end
				# end
			else
				if sp.ext[:options][name][1]
					model.ext[:coverconstraint][node][name] = @constraint(model, 0 >= sum(model.ext[:expansions][past[index]][name] for index in interval))
				else
					model.ext[:coverconstraint][node][name] = @constraint(model, 0 <= sum(model.ext[:expansions][past[index]][name] for index in interval))
				end
				# if typeof(node)==Leaf && sp.ext[:options][name][1]
				# 	@constraint(model, sum(model.ext[:expansions][past][name] for past in history_function(node))<=1)
				# end
			end
		end
	end

	model.ext[:convexcombination] = Dict{AbstractTree,ConstraintRef}()
	for node in keys(sub_problems)
		model.ext[:convexcombination][node] = @constraint(model, 0 == 1)
	end

	if typeof(sideconstraints) <: Function
		map(Main.eval,unpack_expansions(model.ext[:expansions])) #bring expansion variables into global scope
		sideconstraints(model,tree)
		map(Main.eval,clear_expansions(model.ext[:expansions]))
	end

	model
end

function scale_objectives(tree::T where T <: AbstractTree,sub_problems,discount_factor::Float64)
	if discount_factor <= 0.0 || discount_factor > 1.0
		error("discount_factor must be greater than 0.0 and less than or equal to 1.0")
	end

	depth_function=depth(tree)

	for (node,sp) in sub_problems
		if !haskey(sp.ext, :capitalcosts)
			sp.ext[:capitalcosts]=Dict()
		end
		if !haskey(sp.ext, :ongoingcosts)
			sp.ext[:ongoingcosts]=Dict()
		end
		@objective(sp, Min, 0.0)
		set_normalized_coefficient(sp.ext[:objective_con],sp.ext[:objective],1.0/(discount_factor^depth_function(node)))
	end
end

function Base.map(f, hello::JuMP.Containers.DenseAxisArray)
	JuMP.Containers.DenseAxisArray(map(f, hello.data), deepcopy(hello.axes),deepcopy(hello.lookup))
end

function Base.map(f, hello::JuMP.Containers.SparseAxisArray)
	JuMP.Containers.SparseAxisArray(Dict(i => f(hello[i]) for i in keys(hello.data)))
end

function add_variable_as_column(master, info, column)
	holder = JuMP.add_variable(master, JuMP.build_variable(error, info))

	set_normalized_coefficient(master.ext[:convexcombination][column.node], holder, 1.0)

	for var in column.singleVars
		set_normalized_coefficient(master.ext[:coverconstraint][column.node][var], holder, 1.0)
	end

	for (var,array) in column.arrayVars
		for i in array
			set_normalized_coefficient(master.ext[:coverconstraint][column.node][var][i], holder, 1.0)
		end
	end

	for node in collect(column.node)
		if typeof(node)==Leaf
			set_normalized_coefficient(master.ext[:scenprofit_con][node], holder, column.obj)
		end
	end
end

function build_column(master, sub_problem ,node)
	singlevars=Array{Any,1}()
	arrayvars=Dict{Symbol,Array{Any,1}}()
	for (name,variable) in sub_problem.ext[:expansions]
		if typeof(variable) <: AbstractArray
			arrayvars[name]=Array{Int64,1}()
			vals=JuMP.value.(variable)
			for i in eachindex(variable)
				if vals[i] >= 0.99
					push!(arrayvars[name],i)
				end
			end
		else
			if JuMP.value(variable) >= 0.99
				push!(singlevars,name)
			end
		end
	end
	Column(node,singlevars,arrayvars,JuMP.value(sub_problem.ext[:objective]))
end

"""
	solve(judge::JuDGEModel;
	      abstol = 10^-14,
	      reltol = 10^-14,
	      rlx_abstol = 10^-14,
	      rlx_reltol = 10^-14,
	      duration = Inf,
	      iter = 2^63 - 1,
	      inttol = 10^-9,
	      max_no_int = 1000,
	      partial = 100000,
	      warm_starts = false,
	      optimizer_attributes,
	      mp_callback=nothing,
	      allow_frac = :binary_solve,
	      prune = Inf)

Solve a JuDGEModel `judge` without branch and price.

### Required Arguments
`judge` is the JuDGE model that we wish to solve.

### Optional Arguments
`abstol` is the absolute tolerance for the best integer-feasible objective value and the lower bound

`reltol` is the relative tolerance for the best integer-feasible objective value and the lower bound

`rlx_abstol` is the absolute tolerance for the relaxed master objective value and the lower bound

`rlx_reltol` is Set the relative tolerance for the relaxed master objective value and the lower bound

`duration` is the maximum duration

`iter` is the maximum number of iterations

`inttol` is the maximum deviation from 0 or 1 for integer feasible solutions

`max_no_int` is the maximum number of iterations yielding a fractional solution before a MIP solve is
performed on the master

`partial` specifies the number of nodes to solve in each iteration, after all nodes have been solved, a
full pricing iteration is used to compute an updated lower bound

`warm_starts` boolean specifing whether to use warm starts for binary solves of master problem

`optimizer_attributes` can be set to a specific function that dynamically changes optimizer attributes
for the subproblems; this should only be used by people who have examined the source code

`mp_callback` is a user-defined function that specifies termination conditions for MIP solves of the
master problem. See examples/advanced.jl.

### Used by the branch and price algorithm
`allow_frac` indicates whether a fractional solution will be returned; possible values are:
	`:binary_solve` a binary solve of master will be performed (if needed) prior to the solution being returned;
	`:binary_solve_return_relaxation` a binary solve of master will be performed (if needed), updating the upper bound,
	but the master problem relation will be returned;
	`:first_fractional` will return the first fractional master solution found;
	`:no_binary_solve` will simply return the solution to the relaxed master problem when converged.

`prune` is used to stop the algorithm before convergence, if a known upper bound for the problem is specified

### Examples
    JuDGE.solve(jmodel, rlx_abstol=10^-6)
	JuDGE.solve(jmodel, abstol=10^-6)
"""
function solve(judge::JuDGEModel;
   abstol= 10^-14,
   reltol= 10^-14,
   rlx_abstol= 10^-14,
   rlx_reltol= 10^-14,
   duration= Inf,
   iter= 2^63 - 1,
   inttol=10^-9,
   max_no_int=1000,
   partial=100000,
   warm_starts=false,
   optimizer_attributes=nothing,
   mp_callback=nothing,
   allow_frac=:binary_solve,
   prune=Inf,
   verbose=2
   )

	# encode the user convergence test in a ConvergenceState struct
	done = ConvergenceState(0.0, 0.0, 0.0, abstol, reltol, rlx_abstol, rlx_reltol, duration, iter, inttol)
	current = InitialConvergenceState()
	push!(judge.log,current)
	if verbose>0
		print("")
		print("")
		println("Current ObjVal  |   Upper Bound   Lower Bound  |  Absolute Diff   Relative Diff  |  Fractionality  |      Time     Iter")
	end
	# set up times for use in convergence
	initial_time = time()
	obj = Inf
	nodes=collect(judge.tree)
	objduals=Dict{AbstractTree,Float64}()
	redcosts=Dict{AbstractTree,Float64}()
	for node in nodes
		objduals[node]=-1
		redcosts[node]=-1
	end
	no_int_count=0
	if optimizer_attributes!=nothing
		optimizer_attributes(judge,false,true)
	end

	max_char=length(nodes[end].name)+length(string(length(nodes)))
	function get_whitespace(name::String,number::Int64)
		blank="  "
		spaces=max_char-length(name)-length(string(number))
		for i in 1:spaces
			blank*=" "
		end
		blank
	end
	first=-1
	while true
		if first<=0
			nodes2=nodes
		else
			nodes2=AbstractTree[]
			for i in first:min(length(nodes),(first+partial-1))
				push!(nodes2,nodes[i])
			end
		end
		# perform the main iterations
		optimize!(judge.master_problem)
		status=termination_status(judge.master_problem)

		for i in 1:length(nodes2)
			node=nodes2[i]
			updateduals(judge.master_problem, judge.sub_problems[node], node, status, current.iter)
			if verbose==2
				overprint("Solving subproblem for node "*node.name*get_whitespace(node.name,i)*string(i)*"/"*string(length(nodes2)))
			end
			optimize!(judge.sub_problems[node])
			if termination_status(judge.sub_problems[node])!=MathOptInterface.OPTIMAL
				error("Solve for subproblem "*node.name*" exited with status "*string(termination_status(judge.sub_problems[node])))
			end
		end
		if verbose==2
			overprint("")
		end
		frac=NaN
		if status!=MathOptInterface.INFEASIBLE_OR_UNBOUNDED && status!=MathOptInterface.INFEASIBLE && status!=MathOptInterface.DUAL_INFEASIBLE
			if status!=MathOptInterface.OPTIMAL
				@warn("Master problem did not solve to optimality: "*string(status))
			end
			for node in nodes2
				objduals[node]=objective_bound(judge.sub_problems[node])-dual(judge.master_problem.ext[:convexcombination][node])
				redcosts[node]=objective_value(judge.sub_problems[node])-dual(judge.master_problem.ext[:convexcombination][node])
			end
			if first==0
				getlowerbound(judge,objduals)
			end
			frac = absolutefractionality(judge)
			obj = objective_value(judge.master_problem)
			if frac<done.int
				judge.bounds.UB=obj
				if warm_starts
					vars=all_variables(judge.master_problem)
					set_start_value.(vars, JuMP.value.(vars))
				end
				no_int_count=0
				# todo keep track of number of columns added when new UB found
			else
				no_int_count+=1
			end
		elseif judge.bounds.LB>-Inf
			println("\nMaster problem is infeasible or unbounded")
			return
		end

		current = ConvergenceState(obj, judge.bounds.UB, judge.bounds.LB, time() - initial_time, current.iter + 1, frac)
		if verbose>0
			println(current)
		end
		push!(judge.log,current)

		if prune<judge.bounds.LB
			if verbose>0
				println("\nDominated by incumbent.")
			end
			return
		elseif has_converged(done, current)
			solve_master_binary(judge,initial_time,done,allow_frac,warm_starts,nothing,verbose)
			if verbose>0
				println("\nConvergence criteria met.")
			end
			return
		elseif allow_frac==:first_fractional && frac>done.int
			if verbose>0
				println("\nFractional solution found.")
			end
			return
		elseif no_int_count >= max_no_int && current.int > done.int && (current.rlx_abs<done.int_abs || current.rlx_rel<done.int_rel)
			current=solve_master_binary(judge,initial_time,done,:binary_solve_return_relaxation,warm_starts,mp_callback,verbose)
			if has_converged(done, current)
				println("\nConvergence criteria met.")
				return
			end
			no_int_count=0
		end

		num_var=num_variables(judge.master_problem)
		for node in nodes2
			if redcosts[node]<-10^-10
				column = build_column(judge.master_problem, judge.sub_problems[node], node)
				add_variable_as_column(judge.master_problem, UnitIntervalInformation(), column)
				push!(judge.master_problem.ext[:columns],column)
			end
		end

		if optimizer_attributes==nothing
			if first==0 && num_var==num_variables(judge.master_problem)
				solve_master_binary(judge,initial_time,done,allow_frac,warm_starts,nothing,verbose)
				if verbose>0
					println("\nStalled: exiting.")
				end
				return
			end
		elseif optimizer_attributes(judge,num_var==num_variables(judge.master_problem),false)
			if verbose>0
				println("\nStalled: exiting.")
			end
			return
		end
		if first<=0
			first+=1
		else
			first+=partial
		end
		if first+partial>length(nodes)
			first=0
		end
	end
end

function solve_master_binary(judge::JuDGEModel,initial_time::Float64,done::ConvergenceState,allow_frac::Symbol,warm_starts::Bool,mp_callback,verbose)
	set=0
	current=nothing
	if (allow_frac==:binary_solve || allow_frac==:binary_solve_return_relaxation) && judge.log[end].int>done.int
		if verbose==2
			print("Solving master problem as MIP")
		end
		solve_binary(judge,done.int_abs,done.int_rel,warm_starts,mp_callback)
		set=1
		current = ConvergenceState(judge.log[end].obj, judge.bounds.UB, judge.bounds.LB, time()-initial_time, judge.log[end].iter + 1, absolutefractionality(judge))
		if verbose>0
			if verbose==2
				overprint("")
			end
			print(current)
			println("*")
		end
		push!(judge.log,current)
	end
	if allow_frac==:binary_solve_return_relaxation && set==1
		remove_binary(judge)
		optimize!(judge.master_problem)
	end
	current
end

function getlowerbound(judge::JuDGEModel,objduals::Dict{AbstractTree,Float64})
	lb = objective_value(judge.master_problem)
	for (i,objdual) in objduals
		lb+=objdual
	end
	if lb>judge.bounds.LB
		judge.bounds.LB=lb
	end
end

function solve_binary(judge::JuDGEModel,abstol::Float64,reltol::Float64,warm_starts::Bool,mp_callback)
	if typeof(judge.master_solver) <: Tuple
		JuMP.set_optimizer(judge.master_problem,judge.master_solver[2])
	end
	for node in keys(judge.master_problem.ext[:expansions])
		for x in keys(judge.master_problem.ext[:expansions][node])
			var = judge.master_problem.ext[:expansions][node][x]
			if typeof(var) <: AbstractArray
				for key in keys(var)
					set_binary(var[key])
				end
			else
				set_binary(var)
			end
		end
	end
	if typeof(mp_callback) <: Function
		mp_callback(judge,abstol,reltol)
		printright("Initialising")
	end
	optimize!(judge.master_problem)
	obj = objective_value(judge.master_problem)
	if obj<judge.bounds.UB
		judge.bounds.UB=obj
		if warm_starts
			vars=all_variables(judge.master_problem)
			set_start_value.(vars, JuMP.value.(vars))
		end
	end
end

function remove_binary(judge::JuDGEModel)
	for node in keys(judge.master_problem.ext[:expansions])
		for x in keys(judge.master_problem.ext[:expansions][node])
			var = judge.master_problem.ext[:expansions][node][x]
			if typeof(var) <: AbstractArray
				for key in keys(var)
					unset_binary(var[key])
				end
			else
				unset_binary(var)
			end
		end
	end
	if typeof(judge.master_solver) <: Tuple
		JuMP.set_optimizer(judge.master_problem,judge.master_solver[1])
	end
end

function absolutefractionality(jmodel::JuDGEModel;node=jmodel.tree,f=0)
	# this is how you access the value of the binary expansions in the master
	for x in keys(jmodel.master_problem.ext[:expansions][node])
		var = jmodel.master_problem.ext[:expansions][node][x]
		if typeof(var) <: AbstractArray
			for key in keys(var)
				f=max(f,min(JuMP.value(var[key]),1-JuMP.value(var[key])))
			end
		else
			f=max(f,min(JuMP.value(var),1-JuMP.value(var)))
		end
	end

	if typeof(node)==Tree
		for child in node.children
			f=max(f,absolutefractionality(jmodel,node=child,f=f))
		end
	end
	f
end

function updateduals(master, sub_problem, node, status, iter)
	if status!=MathOptInterface.INFEASIBLE_OR_UNBOUNDED && status!=MathOptInterface.INFEASIBLE && status!=MathOptInterface.DUAL_INFEASIBLE
		for (name,var) in sub_problem.ext[:expansions]
			if typeof(var) <: AbstractArray
				for i in eachindex(var)
					set_objective_coefficient(sub_problem, var[i], -dual(master.ext[:coverconstraint][node][name][i]))
				end
			else
				set_objective_coefficient(sub_problem, var, -dual(master.ext[:coverconstraint][node][name]))
			end
		end
		nodes=collect(node)
		total=0.0
		for n in nodes
			if typeof(n)==Leaf
				total-=dual(master.ext[:scenprofit_con][n])
			end
		end
		set_objective_coefficient(sub_problem,sub_problem.ext[:objective],total)
	else
		if iter%2==0
			oc=-10.0^10
		else
			oc=10.0^10
		end
		for (name,var) in sub_problem.ext[:expansions]
			if typeof(var) <: AbstractArray
				for i in eachindex(var)
					set_objective_coefficient(sub_problem, var[i], oc)
				end
			else
				set_objective_coefficient(sub_problem, var, oc)
			end
		end
		set_objective_coefficient(sub_problem,sub_problem.ext[:objective], 1.0)
	end
end

"""
	resolve_subproblems(judge::JuDGEModel)

Once a JuDGE model has converged, it is necessary to re-solve the subproblems to find the optimal decisions within each node.

### Required Arguments
`jmodel` is the JuDGE model that we wish to solve.

### Examples
    resolve_subproblems(judge)
"""
function resolve_subproblems(jmodel::JuDGEModel)
	fix_expansions(jmodel)
	resolve_fixed(jmodel)
end

function fix_expansions(jmodel::JuDGEModel)
	if termination_status(jmodel.master_problem) != MathOptInterface.OPTIMAL && termination_status(jmodel.master_problem) != MathOptInterface.INTERRUPTED
		error("You need to first solve the decomposed model.")
	end

	nodes=collect(jmodel.tree)
	history_fn=history(jmodel.tree)

	for node in nodes
		sp=jmodel.sub_problems[node]
		set_objective_coefficient(sp, sp.ext[:objective], 1.0)

		for (name,var) in jmodel.master_problem.ext[:expansions][node]
			var2=sp[name]
			if typeof(var) <: AbstractArray
				for i in eachindex(var)
					value=0.0
					con_obj=constraint_object(jmodel.master_problem.ext[:coverconstraint][node][name][i])
					for prev in history_fn(node)
						var3=jmodel.master_problem.ext[:expansions][prev][name][i]
						if var3 in keys(con_obj.func.terms)
							value+=JuMP.value(var3)*-con_obj.func.terms[var3]
						end
					end
					if sp.ext[:options][name][1]
						JuMP.set_lower_bound(var2[i],value)
					else
						JuMP.set_upper_bound(var2[i],value)
					end
					set_objective_coefficient(sp, var2[i], 0.0)
				end
			elseif isa(var,VariableRef)
				value=0.0
				con_obj=constraint_object(jmodel.master_problem.ext[:coverconstraint][node][name])
				for prev in history_fn(node)
					var3=jmodel.master_problem.ext[:expansions][prev][name]
					if var3 in keys(con_obj.func.terms)
						value+=JuMP.value(var3)*-con_obj.func.terms[var3]
					end
				end
				if sp.ext[:options][name][1]
					JuMP.set_lower_bound(var2,value)
				else
					JuMP.set_upper_bound(var2,value)
				end
				set_objective_coefficient(sp, var2, 0.0)
			end
		end
	end
end

function resolve_fixed(jmodel::JuDGEModel)
	history_fn=history(jmodel.tree)

	for n in collect(jmodel.tree)
		JuMP.optimize!(jmodel.sub_problems[n])
	end

	scenario_objs=Array{Tuple{Leaf,Float64},1}()
	for (leaf,con) = jmodel.master_problem.ext[:scenprofit_con]
		obj=0.0
		for node in history_fn(leaf)
			obj+=objective_value(jmodel.sub_problems[node])
			for key in keys(jmodel.master_problem.ext[:expansions][node])
				var = jmodel.master_problem.ext[:expansions][node][key]
				if isa(var,VariableRef)
					obj+=JuMP.value(var)*normalized_coefficient(con, var)
				elseif typeof(var) <: AbstractArray
					for v in keys(var)
						obj+=JuMP.value(var[v])*normalized_coefficient(con,var[v])
					end
				end
			end
		end
		push!(scenario_objs,(leaf,obj))
	end

	sort!(scenario_objs,by=i->i[2],rev=true)
	obj=0.0
	beta=jmodel.CVaR[2]
	for scen in scenario_objs
		pr=jmodel.probabilities[scen[1]]
		obj+=scen[2]*pr*(1-jmodel.CVaR[1])
		if pr>beta
			obj+=scen[2]*jmodel.CVaR[1]*beta/jmodel.CVaR[2]
			beta=0
		else
			obj+=scen[2]*jmodel.CVaR[1]*pr/jmodel.CVaR[2]
			beta-=pr
		end
	end
	obj
end

include("model_verification.jl")

# pretty printing
function Base.show(io::IO, ::MIME"text/plain", judge::JuDGEModel)
	print(io, "JuDGE Model with:\n")
	println(io, "  Tree: ", judge.tree)
	print(io, "  Expansion variables: ")
	keys = judge.sub_problems[judge.tree].ext[:expansions]
	for (i,ii) in keys
		print(io, "$(i) ")
	end
end

function Base.show(io::IO, judge::JuDGEModel)
	print(io, "JuDGE Model with:\n")
	println(io, "  Tree: ", judge.tree)
	print(io, "  Expansion variables: ")
	keys = judge.sub_problems[judge.tree].ext[:expansions]
	for (i,ii) in keys
		print(io, "$(i) ")
	end
end

include("output.jl")

export @expansion, @shutdown, @expansionconstraint, @capitalcosts, @ongoingcosts, @sp_objective, JuDGEModel, Leaf, Tree, AbstractTree, narytree, ConditionallyUniformProbabilities, UniformLeafProbabilities, get_node, tree_from_leaves, tree_from_nodes, tree_from_file, DetEqModel, resolve_subproblems

end
