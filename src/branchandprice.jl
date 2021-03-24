function add_branch_constraint(master::JuMP.Model,branch)
	if branch[2]==:eq
		@constraint(master,branch[1]==branch[3])
	elseif branch[2]==:le
		@constraint(master,branch[1]<=branch[3])
	elseif branch[2]==:ge
		@constraint(master,branch[1]>=branch[3])
	end
end

function new_branch_constraint(master::JuMP.Model,branch::Any)
	add_branch_constraint(master,branch)
	push!(master.ext[:branches],branch)
end

function copy_model(jmodel::JuDGEModel, branch::Any, warm_start::Bool)
	newmodel=JuDGEModel(jmodel.tree, jmodel.probabilities, jmodel.sub_problems, jmodel.master_solver, Bounds(Inf,jmodel.bounds.LB), jmodel.discount_factor, jmodel.risk, jmodel.sideconstraints, jmodel.optimizer_settings)
	newmodel.master_problem.ext[:branches]=Array{Any,1}()
	newmodel.master_problem.ext[:sp_branches]=Dict{AbstractTree,Array{Any,1}}()

	for node in collect(jmodel.tree)
		newmodel.master_problem.ext[:sp_branches][node]=Array{Any,1}()
	end

	all_newvar=all_variables(newmodel.master_problem)
	all_var=all_variables(jmodel.master_problem)

	if warm_start
		for i in 1:length(all_newvar)
			MOI.set(newmodel.master_problem, MOI.VariablePrimalStart(), all_newvar[i],
											 MOI.get(jmodel.master_problem, MOI.VariablePrimalStart(),all_var[i]))
		end
	end

	for column in jmodel.master_problem.ext[:columns]
		if typeof(branch)!=Nothing && branch[1]==:subproblem
			if column in branch[5]
				continue
			end
		end
		column2=Column(column.node,column.singleVars,column.arrayVars,column.obj,JuMP.add_variable(newmodel.master_problem, JuMP.build_variable(error, UnitIntervalInformation())),column.solution)
		add_variable_as_column(newmodel.master_problem, column2)
		push!(newmodel.master_problem.ext[:columns],column2)
		push!(all_newvar,column2.var)

		if newmodel.sub_problems[column2.node].ext[:form]==:mixed
			add_mixed_cover(newmodel.master_problem, jmodel.sub_problems[column2.node], column2)
			for i in 1:length(jmodel.sub_problems[column2.node].ext[:discrete])
				push!(all_newvar, jmodel.master_problem.ext[:discrete_var][column2.node][i])
			end
		end

		# MOI.set(newmodel.master_problem, MOI.VariablePrimalStart(), newvar,
		# 								 MOI.get(jmodel.master_problem, MOI.VariablePrimalStart()))
	end

	for b in jmodel.master_problem.ext[:branches]
		if typeof(b[1])==VariableRef
			push!(newmodel.master_problem.ext[:branches],(all_newvar[JuMP.index(b[1]).value],b[2],b[3]))
		elseif typeof(b[1])==GenericAffExpr{Float64,VariableRef}
			exp=0
			for (var,coef) in b[1].terms
				exp+=coef*all_newvar[JuMP.index(var).value]
			end
			push!(newmodel.master_problem.ext[:branches],(exp,b[2],b[3]))
		end
		add_branch_constraint(newmodel.master_problem,newmodel.master_problem.ext[:branches][end])
	end

	if :sp_branches in keys(jmodel.master_problem.ext)
		for (node,sp_branches) in jmodel.master_problem.ext[:sp_branches]
			for b in sp_branches
				push!(newmodel.master_problem.ext[:sp_branches][node],b)
			end
		end
	end

	if typeof(branch)!=Nothing
		if branch[1]==:master
			if typeof(branch[2])==VariableRef
				branch=(all_newvar[JuMP.index(branch[2]).value],branch[3],branch[4])
			elseif typeof(branch[2])==GenericAffExpr{Float64,VariableRef}
				exp=0
				for (var,coef) in branch[2].terms
					exp+=coef*all_newvar[JuMP.index(var).value]
				end
				branch=(exp,branch[3],branch[4])
			end
			new_branch_constraint(newmodel.master_problem,branch)
		elseif branch[1]==:subproblem
			if typeof(branch[2])!=VariableRef
				@error("subproblem branches must be on a single variable")
			end
			newbranch=(branch[2],branch[3],branch[4])
			push!(newmodel.master_problem.ext[:sp_branches][branch[6]],newbranch)
		end
	end

	newmodel
end

function copy_column(master, sub_problem ,node, master_var)
	singlevars=Array{Any,1}()
	arrayvars=Dict{Any,Array{Any,1}}()

	for (name,variable) in sub_problem.ext[:expansions]
		if typeof(variable) <: AbstractArray
			arrayvars[name]=Array{Int64,1}()
			for i in eachindex(variable)
				if normalized_coefficient(master.ext[:coverconstraint][node][name][i],master_var) == 1.0
					push!(arrayvars[name],i)
				end
			end
		else
			if normalized_coefficient(master.ext[:coverconstraint][node][name],master_var) == 1.0
				push!(singlevars,name)
			end
		end
	end
	Column(node,singlevars,arrayvars,objcoef(master_var))
end

function variable_branch_most_frac(master, tree, expansions, inttol)
	branches=Array{Any,1}()

	branch=nothing
	maxfrac=inttol
	for node in collect(tree)
		for x in keys(expansions[node])
			var = expansions[node][x]
			if typeof(var) <: AbstractArray
				for key in keys(var)
					fractionality=min(JuMP.value(var[key]),1-JuMP.value(var[key]))
					if fractionality>maxfrac
						branch=@expression(master,var[key])
						maxfrac=fractionality
					end
				end
			else
				fractionality=min(JuMP.value(var),1-JuMP.value(var))
				if fractionality>maxfrac
					branch=@expression(master,var)
					maxfrac=fractionality
				end
			end
		end
	end
	if branch!=nothing
		push!(branches,(branch,:eq,1))
		push!(branches,(branch,:eq,0))
	end
	branches
end

"""
	variable_branch(master, subproblems, tree, expansions, inttol)

This is an in-built function that is called during branch-and-price to perform a branch.
Users can define their own functions that follow this format to create new branching strategies.

### Required Arguments
`master` is the master problem of the JuDGE model

`subproblems` is a dictionary of the subproblems of the JuDGE model

`tree` is the tree of the JuDGE model

`expansions` is a dictionary of capacity expansion/shutdown variables, indexed by node.

`inttol` is the maximum permitted deviation from 0 or 1 for a value to still be considered binary.
"""
function variable_branch(master, subproblems, tree, expansions, inttol)
	branches=Array{Any,1}()

	for node in collect(tree,order=:breadth)
		if subproblems[node].ext[:form]==:mixed
			cutoff=inttol
			index=0

			for i in 1:length(subproblems[node].ext[:discrete_branch])
				val=JuMP.value(master.ext[:discrete_var][node][i])
				fractionality=min(val-floor(val),ceil(val)-val)
				if fractionality>inttol
					cutoff=val
					index=i
				end
			end

			if index!=0
				group1=[]
				group2=[]

				for col in master.ext[:columns]
					if col.node==node
						val=col.solution[index]

						if val<cutoff
							push!(group2,col)
						else
							push!(group1,col)
						end
					end
				end
				branch=@expression(subproblems[node],subproblems[node].ext[:discrete_branch][index])

				push!(branches,(:subproblem,branch,:le,ceil(cutoff)-1,group1,node))
				push!(branches,(:subproblem,branch,:ge,ceil(cutoff),group2,node))
				return branches
			end
		end
		for x in keys(expansions[node])
			if master.ext[:options][x][4]!=:Con
				var = expansions[node][x]
				if typeof(var) <: AbstractArray
					if typeof(var) <: JuMP.Containers.SparseAxisArray
						var=var.data
					end
					for key in keys(var)
						val=JuMP.value(var[key])
						fractionality=min(val-floor(val),ceil(val)-val)
						if fractionality>inttol
							branch=@expression(master,var[key])
							push!(branches,(:master,branch,:ge,ceil(val)))
							push!(branches,(:master,branch,:le,floor(val)))
							return branches
						end
					end
				else
					val=JuMP.value(var)
					fractionality=min(val-floor(val),ceil(val)-val)
					if fractionality>inttol
						branch=@expression(master,var)
						push!(branches,(:master,branch,:ge,ceil(val)))
						push!(branches,(:master,branch,:le,floor(val)))
						return branches
					end
				end
			end
		end
	end
	branches
end

function perform_branch(jmodel::JuDGEModel,branches::Array{Any,1}, warm_starts::Bool)
	newmodels=Array{JuDGEModel,1}()

	for i in 1:length(branches)
		push!(newmodels,copy_model(jmodel,branches[i],warm_starts))
	end

	#new_branch_constraint(jmodel.master_problem,branches[1])
	#jmodel.bounds.UB=Inf
	#insert!(newmodels,1,jmodel)

	newmodels
end

"""
	branch_and_price(judge::JuDGEModel;
		  branch_method=JuDGE.variable_branch,
		  search=:depth_first_resurface,
	      abstol = 10^-14,
	      reltol = 10^-14,
	      rlx_abstol = 10^-14,
	      rlx_reltol = 10^-14,
	      duration = Inf,
	      iter = 2^63 - 1,
	      inttol = 10^-7,
	      max_no_int = 1000,
	      warm_starts = false,
	      allow_frac=:binary_solve_return_relaxation,
		  optimizer_attributes,
	      mp_callback=nothing)

Solve a JuDGEModel `judge` without branch and price.

### Required Arguments
`judge` is the JuDGE model that we wish to solve.

### Optional Arguments
`branch_method` specifies the way that constrants are added to create new nodes

`search` specifies the order in which nodes are solved in the (brnach-and-price) tree

`abstol` is the absolute tolerance for the best integer-feasible objective value and the lower bound

`reltol` is the relative tolerance for the best integer-feasible objective value and the lower bound

`rlx_abstol` is the absolute tolerance for the relaxed master objective value and the lower bound

`rlx_reltol` is Set the relative tolerance for the relaxed master objective value and the lower bound

`duration` is the maximum duration

`iter` is the maximum number of iterations

`inttol` is the maximum deviation from 0 or 1 for integer feasible solutions

`max_no_int` is the maximum number of iterations yielding a fractional solution before a MIP solve is
performed on the master

`warm_starts` boolean specifing whether to use warm starts for binary solves of master problem

`allow_frac` indicates whether a fractional solution will be returned by `JuDGE.solve()`; possible values are:
	`:binary_solve_return_relaxation` a binary solve of master will be performed (if needed), updating the upper bound,
	but the master problem relation will be returned by `JuDGE.solve()`;
	`:first_fractional` `JuDGE.solve()` will return the first fractional master solution found;
	`:no_binary_solve` `JuDGE.solve()` will simply return the solution to the relaxed master problem when converged.

`partial` specifies the number of nodes to solve in each iteration, after all nodes have been solved, a
full pricing iteration is used to compute an updated lower bound

`optimizer_attributes` can be set to a specific function that dynamically changes optimizer attributes
for the subproblems; this should only be used by people who have examined the source code

`mp_callback` is a user-defined function that specifies termination conditions for MIP solves of the
master problem. See examples/advanced.jl.

### Examples
	JuDGE.branch_and_price(jmodel, abstol=10^-6)
	JuDGE.branch_and_price(jmodel, branch_method=JuDGE.variable_branch, search=:lowestLB)
"""
function branch_and_price(models::Union{JuDGEModel,Array{JuDGEModel,1}};branch_method=JuDGE.variable_branch,search=:lowestLB,
   abstol=10^-14,
   reltol=10^-14,
   rlx_abstol=10^-14,
   rlx_reltol=10^-14,
   duration=Inf,
   iter=2^63-1,
   inttol=10^-7,
   max_no_int=1000,
   warm_starts=false,
   allow_frac=:binary_solve_return_relaxation,
   blocks=nothing,
   verbose=2,
   optimizer_attributes=nothing,
   mp_callback=nothing)
#	if judge.master_problem.ext[:types]!=:binary
#		@warn("The branch and price implementation is not guaranteed to converge with continuous expansions.")
#	end
	initial_time = time()
	#log=[]

	if typeof(rlx_abstol) <: Function
		rlx_abstol_func=true
	elseif typeof(rlx_abstol)==Float64
		rat=rlx_abstol
		rlx_abstol_func=false
	end
	if typeof(rlx_reltol) <: Function
		rlx_reltol_func=true
	elseif typeof(rlx_reltol)==Float64
		rrt=rlx_reltol
		rlx_reltol_func=false
	end

	if typeof(models)==JuDGEModel
		models = [models]
	end
	#push!(models,copy_model(judge,nothing))
	#push!(models,judge)
	models[1].master_problem.ext[:branches]=Array{Any,1}()

	UB=Inf
	LB=Inf
	otherLB=Inf
	i=1
	best=models[1]
	while true
		model=models[i]

		N=length(models)
		while model.bounds.LB>UB
			if verbose>0
				print("\n")
			end
			println("Model "*string(i)*" dominated. UB: "*string(model.bounds.UB)*", LB:"*string(model.bounds.LB))
			if i==N
				break
			end
			i+=1
			model=models[i]
		end

		bestLB=0.0
		LB=Inf
		for j in i:N
			if models[j].bounds.LB<LB
				LB=models[j].bounds.LB
				bestLB=j
			end
		end
		if otherLB<LB
			LB=otherLB
		end

		if search==:lowestLB
			model=models[bestLB]
			deleteat!(models,bestLB)
			insert!(models,i,model)
		end

		if verbose>0
			print("\n")
		end
		println("Model "*string(i)*" of "*string(N)*". UB: "*string(UB)*", LB:"*string(LB)*", Time: "*string(Int(floor((time()-initial_time)*1000+0.5))/1000)*"s")
		#push!(log,(UB,LB,(time()-initial_time)))
		if rlx_reltol_func
			rrt=rlx_reltol(i,N)
		end

		if rlx_abstol_func
			rat=rlx_abstol(i,N)
		end

		solve(model,abstol=abstol,reltol=reltol,rlx_abstol=rat,rlx_reltol=rrt,duration=duration,iter=iter,inttol=inttol,warm_starts=warm_starts,allow_frac=allow_frac,prune=UB,optimizer_attributes=optimizer_attributes,mp_callback=mp_callback,max_no_int=max_no_int,blocks=blocks,verbose=verbose)

		if model.bounds.UB<UB
			UB=model.bounds.UB
			best=copy_model(model,nothing,warm_starts)
		end

		bestLB=0
		LB=otherLB
		for j in i:N
			if models[j].bounds.LB<LB
				LB=models[j].bounds.LB
				bestLB=j
			end
		end

		status=termination_status(model.master_problem)
		if model.bounds.LB<=UB && (LB+abstol<UB && UB-LB>reltol*abs(LB)) && (status!=MathOptInterface.INFEASIBLE_OR_UNBOUNDED && status!=MathOptInterface.INFEASIBLE && status!=MathOptInterface.DUAL_INFEASIBLE && status!=MathOptInterface.NUMERICAL_ERROR)
			if verbose>0
				println("\nAttempting to branch.")
			end
			branches=branch_method(model.master_problem, model.sub_problems, model.tree, model.master_problem.ext[:expansions], inttol)
			if length(branches)>0
				if verbose>0
					println("Adding "*string(length(branches))*" new nodes to B&P tree.")
				end
			    newmodels=perform_branch(model,branches,warm_starts)
				i+=1
				if search==:depth_first_dive
					for j in 1:length(newmodels)
						insert!(models,i,newmodels[length(newmodels)+1-j])
					end
				elseif search==:breadth_first || search==:lowestLB
					append!(models,newmodels)
				elseif search==:depth_first_resurface
					insert!(models,i,newmodels[1])
					deleteat!(newmodels,1)
					append!(models,newmodels)
				end
			elseif i==length(models)
				break
			else
				if model.bounds.LB<otherLB
					otherLB=model.bounds.LB
				end
				i+=1
			end
		elseif i==length(models) || UB-LB<=abstol || UB-LB<=reltol*abs(UB)
			break
		else
			if model.bounds.LB<otherLB
				otherLB=model.bounds.LB
			end
			i+=1
		end
	end

	if allow_frac!=:no_binary_solve
		if verbose==2
			print("Performing final MIP solve")
		end
		solve_binary(best,abstol,reltol,warm_starts,nothing)
		if verbose==2
			overprint("")
		end
	else
		optimize!(best.master_problem)
	end
	UB=best.bounds.UB
	#push!(log,(UB,LB,(time()-initial_time)))
	println("\nObjective value of best integer-feasible solution: "*string(UB))
	println("Objective value of lower bound: "*string(min(LB,UB)))
	println("Solve time: "*string(Int(floor((time()-initial_time)*1000+0.5))/1000)*"s")
	best#, log
end
