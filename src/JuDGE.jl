module JuDGE

using MathOptInterface
using Printf
using Distributed
@everywhere using JuMP

include("tree.jl")
include("macros.jl")
include("convergence.jl")
include("risk.jl")
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
	var::VariableRef
	solution::Array{Float64,1}
end

const RiskNeutral=nothing

struct JuDGEModel
	tree::AbstractTree
	master_problem::JuMP.Model
	sub_problems::Dict{AbstractTree,JuMP.Model}
	bounds::Bounds
	discount_factor::Float64
	master_solver
	probabilities
	risk::Union{Nothing,Risk,Array{Risk,1}}
	sideconstraints
	log::Array{ConvergenceState,1}
	optimizer_settings::Array{Int64,1}
#	variable_bounds::Dict{VariableRef,Tuple{Float64,Float64}}
end

include("utilities.jl")

function JuDGEModel(tree::T where T <: AbstractTree, probabilities::Dict{AbstractTree,Float64}, sub_problems::Dict{AbstractTree,JuMP.Model}, solver, bounds::Bounds, discount_factor::Float64, risk::Union{Nothing,Risk,Array{Risk,1}}, sideconstraints, optimizer_settings::Array{Int64,1})
	master_problem = build_master(sub_problems, tree, probabilities, solver, discount_factor, risk, sideconstraints)
	JuDGEModel(tree,master_problem,sub_problems,bounds,discount_factor,solver,probabilities,risk,sideconstraints,Array{ConvergenceState,1}(),optimizer_settings)
end

"""
	JuDGEModel(tree::AbstractTree,
               probabilities,
               sub_problem_builder::Function,
               solver;
               discount_factor=1.0,
               risk=RiskNeutral,
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

`risk` can be tuple with the two CVaR parameters: (λ, α), but also can be specified as a vector of such tuples.

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
                                    discount_factor=0.9, risk=(0.5,0.1)))
"""
function JuDGEModel(tree::T where T <: AbstractTree, probabilities, sub_problem_builder::Function, solver; discount_factor=1.0, risk::Union{Nothing,Risk,Array{Risk,1}}=nothing, sideconstraints=nothing, parallel=false, sp_solver=nothing, check=true, perfect_foresight=false)
	println("")
	if !perfect_foresight
		println("Establishing JuDGE model for tree: " * string(tree))
	else
		println("Establishing perfect foresight models for tree: " * string(tree))
	end

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

	if !perfect_foresight
		print("Building master problem...")
		master_problem = build_master(sub_problems, tree, probabilities, solver, discount_factor, risk, sideconstraints)
		println("Complete")
		return JuDGEModel(tree,master_problem,sub_problems,Bounds(Inf,-Inf),discount_factor,solver,probabilities,risk,sideconstraints,Array{ConvergenceState,1}(),Int64[])
	else
		scenarios=Dict{AbstractTree,JuDGEModel}()
		pr=Dict{AbstractTree,Float64}()
		print("Building master problems...")
		scen_trees=get_scenarios(tree)
		for t in scen_trees
			leaf=getID(get_leafnodes(t)[1])
			sps = Dict(i => sub_problems[getID(i)] for i in collect(t))
			probs = Dict(i => 1.0 for i in collect(t))
			master_problem = build_master(sps, t, probs, solver, discount_factor, risk, sideconstraints)
			scenarios[leaf]=JuDGEModel(t,master_problem,sps,Bounds(Inf,-Inf),discount_factor,solver,probabilities,risk,sideconstraints,Array{ConvergenceState,1}(),Int64[])
			pr[leaf]=probabilities[leaf]
		end
		println("Complete")
		return scenarios, pr
	end
end

include("branchandprice.jl")
include("master.jl")

function scale_objectives(tree::T where T <: AbstractTree,sub_problems,discount_factor::Float64)
	if discount_factor <= 0.0 || discount_factor > 1.0
		error("discount_factor must be greater than 0.0 and less than or equal to 1.0")
	end

	for (node,sp) in sub_problems
		if !haskey(sp.ext, :capitalcosts)
			sp.ext[:capitalcosts]=Dict()
			sp.ext[:capitalcosts][:constant]=0
		end
		if !haskey(sp.ext, :ongoingcosts)
			sp.ext[:ongoingcosts]=Dict()
			sp.ext[:ongoingcosts][:constant]=0
		end

		sp.ext[:objective]=@variable(sp, obj)
		sp.ext[:objective_con]=@constraint(sp, sp.ext[:objective]-objective_function(sp) == 0)

		@objective(sp, Min, 0.0)
		set_normalized_coefficient(sp.ext[:objective_con],sp.ext[:objective],1.0/(discount_factor^depth(node)))
	end
end

function Base.map(f, hello::JuMP.Containers.DenseAxisArray)
	JuMP.Containers.DenseAxisArray(map(f, hello.data), deepcopy(hello.axes),deepcopy(hello.lookup))
end

function Base.map(f, hello::JuMP.Containers.SparseAxisArray)
	JuMP.Containers.SparseAxisArray(Dict(i => f(hello[i]) for i in keys(hello.data)))
end

function add_variable_as_column(master, column)
	for constr in master.ext[:convexcombination][column.node]
		set_normalized_coefficient(constr, column.var, 1.0)
	end

	for (var,val) in column.singleVars
		set_normalized_coefficient(master.ext[:coverconstraint][column.node][var], column.var, val)
	end

	for (var,array) in column.arrayVars
		for (i,val) in array
			set_normalized_coefficient(master.ext[:coverconstraint][column.node][var][i], column.var, val)
		end
	end

	for node in collect(column.node)
		if typeof(node)==Leaf
			set_normalized_coefficient(master.ext[:scenprofit_con][node], column.var, column.obj)
		end
	end
end

function add_mixed_cover(master, sp, column)
	if !(column.node in keys(master.ext[:discrete_con]))
		master.ext[:discrete_var][column.node]=Dict{Int64,VariableRef}()
		master.ext[:discrete_con][column.node]=Dict{Int64,ConstraintRef}()
		for i in 1:length(sp.ext[:discrete])
			master.ext[:discrete_var][column.node][i]=@variable(master)
			master.ext[:discrete_con][column.node][i]=@constraint(master,0.0==master.ext[:discrete_var][column.node][i])
		end
	end

	for i in 1:length(sp.ext[:discrete])
		set_normalized_coefficient(master.ext[:discrete_con][column.node][i],column.var,column.solution[i])
	end
end

#	h=hash(column.solution[1])
#	if !(h in keys(master.ext[:mixed_cover_con][column.node]))
#		master.ext[:mixed_cover_var][column.node][h]=@variable(master)

#		master.ext[:mixed_cover_con][column.node][h]=@constraint(master,0.0==master.ext[:mixed_cover_var][column.node][h])
#
#	master.ext[:mixed_solution][column.node][h]=column.solution[1]									#
#	end
#
#	set_normalized_coefficient(master.ext[:mixed_cover_con][column.node][h], column.var, 1.0)
#
#	master.ext[:mixed_cover_var][column.node][h]
#end

function build_column(master, sub_problem ,node, sol)
	singlevars=Array{Any,1}()
	arrayvars=Dict{Symbol,Array{Any,1}}()
	for (name,variable) in sub_problem.ext[:expansions]
		if typeof(variable) <: AbstractArray
			arrayvars[name]=Array{Int64,1}()
			vals=JuMP.value.(variable)
			for i in eachindex(variable)
				if sub_problem.ext[:options][name][4]!=:Con
					rval=round(vals[i])
					if rval!=0.0
						push!(arrayvars[name],(i,rval))
					end
				else
					push!(arrayvars[name],(i,vals[i]))
				end
			end
		else
			if sub_problem.ext[:options][name][4]!=:Con
				rval=round(JuMP.value(variable))
				if rval!=0.0
					push!(singlevars,(name,rval))
				end
			else
				push!(singlevars,(name,JuMP.value(variable)))
			end
		end
	end
	Column(node,singlevars,arrayvars,JuMP.value(sub_problem.ext[:objective]),JuMP.add_variable(master, JuMP.build_variable(error, UnitIntervalInformation())),sol)
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

`blocks` specifies the groups of nodes to solve in each iteration, after all nodes have been solved, a
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
	`:no_binary_solve` will simply return the solution to the relaxed master problem when terminated.

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
   blocks=nothing,
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
	if blocks==nothing
		blocks=[nodes]
	end

	block=-1
	# remove_binary(judge)
	optimize!(judge.master_problem)
	status=termination_status(judge.master_problem)
	if status==MOI.NUMERICAL_ERROR
		println("\nMaster problem returned a MOI.NUMERICAL_ERROR")
		return
	end

	while true
		if block<=0
			nodes2 = nodes
		else
			nodes2 = blocks[block]
		end
		b_con=Array{Any,1}()
		# perform the main iterations
		for i in 1:length(nodes2)
			node=nodes2[i]
			sp=judge.sub_problems[node]
			updateduals(judge.master_problem, sp, node, status, current.iter)
			if verbose==2
				overprint("Solving subproblem for node "*node.name*get_whitespace(node.name,i)*string(i)*"/"*string(length(nodes2)))
			end

			if :sp_branches in keys(judge.master_problem.ext)
				if node in keys(judge.master_problem.ext[:sp_branches])
					sp_branches = judge.master_problem.ext[:sp_branches][node]

					for b in sp_branches
						con=nothing
						if b[2]==:eq
							con=@constraint(sp,b[1]==b[3])
						elseif b[2]==:le
							con=@constraint(sp,b[1]<=b[3])
						elseif b[2]==:ge
							con=@constraint(sp,b[1]>=b[3])
						end
						push!(b_con,(sp,con))
					end
				end
			end

			optimize!(sp)

			if termination_status(sp)!=MathOptInterface.OPTIMAL && termination_status(sp)!=MathOptInterface.INTERRUPTED && termination_status(sp)!=MathOptInterface.TIME_LIMIT
				@warn("Solve for subproblem "*node.name*" exited with status "*string(termination_status(sp)))
				println("Subproblem is infeasible or unbounded")
				for (sp,con) in b_con
					delete(sp,con)
				end
				return :sp_infeasible
			end
			if warm_starts
				vars=all_variables(sp)
				set_start_value.(vars, JuMP.value.(vars))
			end
		end
		if verbose==2
			overprint("")
		end
		frac=NaN
		if status!=MathOptInterface.INFEASIBLE_OR_UNBOUNDED && status!=MathOptInterface.INFEASIBLE && status!=MathOptInterface.DUAL_INFEASIBLE
			for node in nodes2
				objduals[node]=objective_bound(judge.sub_problems[node])
				redcosts[node]=objective_value(judge.sub_problems[node])
			end
			if status!=MathOptInterface.OPTIMAL
				@warn("Master problem did not solve to optimality: "*string(status))
			elseif block==0
				getlowerbound(judge,objduals)
			end
		elseif judge.bounds.LB>-Inf && length(b_con)==0
			println("\nMaster problem is infeasible or unbounded")
			for (sp,con) in b_con
				delete(sp,con)
			end
			return
		end

		# if warm_starts && (status!=MathOptInterface.INFEASIBLE_OR_UNBOUNDED && status!=MathOptInterface.INFEASIBLE && status!=MathOptInterface.DUAL_INFEASIBLE)
		# 	vars=all_variables(judge.master_problem)
		# 	set_start_value.(vars, JuMP.value.(vars))
		# end

		num_var=num_variables(judge.master_problem)
		for node in nodes2
			if redcosts[node]<-10^-10 || status==MathOptInterface.INFEASIBLE_OR_UNBOUNDED || status==MathOptInterface.INFEASIBLE || status==MathOptInterface.DUAL_INFEASIBLE
				if judge.sub_problems[node].ext[:form]==:mixed
					sol=JuMP.value.(judge.sub_problems[node].ext[:discrete])
				else
					sol=[0.0]
				end
				column = build_column(judge.master_problem, judge.sub_problems[node], node, sol)
				if warm_starts
					set_start_value(column.var,0.0)
				end
				add_variable_as_column(judge.master_problem, column)
				push!(judge.master_problem.ext[:columns],column)
				if judge.sub_problems[node].ext[:form]==:mixed
					add_mixed_cover(judge.master_problem, judge.sub_problems[node], column)
				end
			end
		end

		for (sp,con) in b_con
			delete(sp,con)
		end

		optimize!(judge.master_problem)
		status=termination_status(judge.master_problem)
		if status==MOI.NUMERICAL_ERROR
			println("\nMaster problem returned a MOI.NUMERICAL_ERROR")
			return
		elseif status!=MathOptInterface.INFEASIBLE_OR_UNBOUNDED && status!=MathOptInterface.INFEASIBLE && status!=MathOptInterface.DUAL_INFEASIBLE
			frac = absolutefractionality(judge)
			obj = objective_value(judge.master_problem)
			if frac<done.int
				judge.bounds.UB=obj
				no_int_count=0
			else
				no_int_count+=1
			end
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
		elseif (no_int_count >= max_no_int && current.int > done.int && (current.rlx_abs<done.int_abs || current.rlx_rel<done.int_rel)) || (max_no_int<0 && no_int_count>-max_no_int)
			current=solve_master_binary(judge,initial_time,done,allow_frac,warm_starts,mp_callback,verbose)
			if has_converged(done, current)
				if verbose>0
					println("\nConvergence criteria met.")
				end
				return
			elseif allow_frac==:binary_solve
				remove_binary(judge)
				optimize!(judge.master_problem)
			end
			no_int_count=0
		end

		if optimizer_attributes==nothing
			if block==0 && num_var==num_variables(judge.master_problem)
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
		if length(blocks)==1
			block=0
		else
			block = (block+1) % (length(blocks)+1)
		end
	end
end

function getlowerbound(judge::JuDGEModel,objduals::Dict{AbstractTree,Float64})
	lb = objective_value(judge.master_problem)
	for (i,objdual) in objduals
		if :le in keys(i.ext)
			lb+=objdual*i.ext[:le]
		else
			lb+=objdual
		end
	end
	if lb>judge.bounds.LB
		judge.bounds.LB=lb
	end
end

function absolutefractionality(jmodel::JuDGEModel;node=jmodel.tree,f=0)
	# this is how you access the value of the binary expansions in the master
	for x in keys(jmodel.master_problem.ext[:expansions][node])
		if jmodel.sub_problems[jmodel.tree].ext[:options][x][4]!=:Con
			var = jmodel.master_problem.ext[:expansions][node][x]
			if typeof(var) <: AbstractArray
				for key in eachindex(var)
					val=JuMP.value(var[key])
					f=max(f,min(val-floor(val),ceil(val)-val))
				end
			else
				val=JuMP.value(var)
				f=max(f,min(val-floor(val),ceil(val)-val))
			end
		end
	end
	if jmodel.sub_problems[node].ext[:form]==:mixed
		for (i,var) in jmodel.master_problem.ext[:discrete_var][node]
			val=JuMP.value(var)
			f=max(f,min(val-floor(val),ceil(val)-val))
		end
	end

	if jmodel.master_problem.ext[:binarycolumns]
		for col in jmodel.master_problem.ext[:columns]
			val=JuMP.value(col.var)
			f=max(f,min(val-floor(val),ceil(val)-val))
		end
	else
		if typeof(node)==Tree
			for child in node.children
				f=max(f,absolutefractionality(jmodel,node=child,f=f))
			end
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

		cc_sum=0.0
		for constr in master.ext[:convexcombination][node]
			set=typeof(constraint_object(constr).set)
			if set <: MathOptInterface.LessThan
				cc_sum+=dual(constr)#*node.ext[:le]
			elseif set <: MathOptInterface.GreaterThan
				cc_sum-=dual(constr)#*node.ext[:ge]
			elseif set <: MathOptInterface.EqualTo
				# mult=1.0
				# if :eq in keys(node.ext)
				# 	mult=node.ext[:eq]
				# end
				cc_sum+=dual(constr)#*mult
			end
		end

		set_objective_function(sub_problem,
					objective_function(sub_problem)-objective_function(sub_problem).constant
					                               -cc_sum)
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

function get_objval(jmodel::JuDGEModel; risk=jmodel.risk)
	scenario_objs=Dict{Leaf,Float64}()

	for (leaf, var) in jmodel.master_problem.ext[:scenprofit_var]
		scenario_objs[leaf]=JuMP.value(var)
	end

	compute_objval(scenario_objs, jmodel.probabilities, risk)
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
	if termination_status(jmodel.master_problem) != MathOptInterface.OPTIMAL && termination_status(jmodel.master_problem) != MathOptInterface.INTERRUPTED && termination_status(jmodel.master_problem) != MathOptInterface.LOCALLY_SOLVED  && termination_status(jmodel.master_problem) != MathOptInterface.INTERRUPTED
		error("You need to first solve the decomposed model.")
	end

	for node in collect(jmodel.tree)
		sp=jmodel.sub_problems[node]
		set_objective_function(sp, objective_function(sp)-objective_function(sp).constant)
		set_objective_coefficient(sp, sp.ext[:objective], 1.0)
		for (name,var) in jmodel.master_problem.ext[:expansions][node]
			var2=sp.ext[:expansions][name]
			if typeof(var) <: AbstractArray
				for i in eachindex(var)
					value=0.0
					if sp.ext[:options][name][1]==:state
						prev = node.parent
						if prev==nothing
							var3=jmodel.master_problem.ext[:expansions][node][name][i]
							value=JuMP.value(var3)-sp.ext[:options][name][7]
						else
							var3=jmodel.master_problem.ext[:expansions][node][name][i]
							var4=jmodel.master_problem.ext[:expansions][prev][name][i]
							value=JuMP.value(var3)-JuMP.value(var4)
						end
					else
						con_obj=constraint_object(jmodel.master_problem.ext[:coverconstraint][node][name][i])
						for prev in history(node)
							var3=jmodel.master_problem.ext[:expansions][prev][name][i]
							if var3 in keys(con_obj.func.terms)
								value+=JuMP.value(var3)*-con_obj.func.terms[var3]
							end
						end
					end
					if sp.ext[:options][name][1]==:shutdown
						JuMP.set_lower_bound(var2[i],value)
					elseif sp.ext[:options][name][1]==:expansion
						JuMP.set_upper_bound(var2[i],value)
					elseif sp.ext[:options][name][1]==:enforced || sp.ext[:options][name][1]==:state
						JuMP.fix(var2[i],value,force=true)
					end
					set_objective_coefficient(sp, var2[i], 0.0)
				end
			elseif isa(var,VariableRef)
				value=0.0
				if sp.ext[:options][name][1]==:state
					prev = node.parent
					if prev==nothing
						var3=jmodel.master_problem.ext[:expansions][node][name]
						value=JuMP.value(var3)-sp.ext[:options][name][7]
					else
						var3=jmodel.master_problem.ext[:expansions][node][name]
						var4=jmodel.master_problem.ext[:expansions][prev][name]
						value=JuMP.value(var3)-JuMP.value(var4)
					end
				else
					con_obj=constraint_object(jmodel.master_problem.ext[:coverconstraint][node][name])
					for prev in history(node)
						var3=jmodel.master_problem.ext[:expansions][prev][name]
						if var3 in keys(con_obj.func.terms)
							value+=JuMP.value(var3)*-con_obj.func.terms[var3]
						end
					end
				end
				if sp.ext[:options][name][1]==:shutdown
					JuMP.set_lower_bound(var2,value)
				elseif sp.ext[:options][name][1]==:expansion
					JuMP.set_upper_bound(var2,value)
				elseif sp.ext[:options][name][1]==:enforced || sp.ext[:options][name][1]==:state
					JuMP.fix(var2,value,force=true)
				end
				set_objective_coefficient(sp, var2, 0.0)
			end
		end
	end
end

function resolve_fixed(jmodel::JuDGEModel)
	for n in collect(jmodel.tree)
		JuMP.optimize!(jmodel.sub_problems[n])
	end

	scenario_objs=Dict{Leaf,Float64}()
	for (leaf,con) = jmodel.master_problem.ext[:scenprofit_con]
		obj=0.0
		for node in history(leaf)
			obj+=objective_value(jmodel.sub_problems[node])
			for key in keys(jmodel.master_problem.ext[:expansions][node])
				var = jmodel.master_problem.ext[:expansions][node][key]
				if isa(var,VariableRef)
					obj+=JuMP.value(var)*normalized_coefficient(con, var)
				elseif typeof(var) <: AbstractArray
					for v in eachindex(var)
						obj+=JuMP.value(var[v])*normalized_coefficient(con,var[v])
					end
				end
			end
		end
		scenario_objs[leaf]=obj;
	end

	compute_objval(scenario_objs, jmodel.probabilities, jmodel.risk)
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

export @expansion, @shutdown, @enforced, @state, @expansionconstraint, @capitalcosts, @ongoingcosts, JuDGEModel, Risk, Leaf, Tree, AbstractTree, narytree, ConditionallyUniformProbabilities, UniformLeafProbabilities, get_node, tree_from_leaves, tree_from_nodes, tree_from_file, DetEqModel, resolve_subproblems

end
