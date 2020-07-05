function add_branch_constraint(master::JuMP.Model,branch)
	# if typeof(branch[1])==VariableRef
	# 	if branch[2]==:eq
	# 		JuMP.fix(branch[1],branch[3],force=true)
	# 	elseif branch[2]==:le
	# 		JuMP.set_upper_bound(branch[1],branch[3])
	# 	elseif branch[2]==:ge
	# 		JuMP.set_lower_bound(branch[1],branch[3])
	# 	end
	# else
		if branch[2]==:eq
			@constraint(master,branch[1]==branch[3])
		elseif branch[2]==:le
			@constraint(master,branch[1]<=branch[3])
		elseif branch[2]==:ge
			@constraint(master,branch[1]>=branch[3])
		end
	# end
end

function new_branch_constraint(master::JuMP.Model,branch::Any)
	add_branch_constraint(master,branch)
	push!(master.ext[:branches],branch)
end

function copy_model(jmodel::JuDGEModel, branch::Any)
	newmodel=JuDGEModel(jmodel.tree, ConditionallyUniformProbabilities, jmodel.sub_problems, jmodel.master_solver, Bounds(Inf,jmodel.bounds.LB), jmodel.discount_factor)

	newmodel.master_problem.ext[:branches]=Array{Any,1}()

	all_newvar=all_variables(newmodel.master_problem)
	all_var=all_variables(jmodel.master_problem)
	num=length(all_newvar)

	for branch in jmodel.master_problem.ext[:branches]
		if typeof(branch[1])==VariableRef
			push!(newmodel.master_problem.ext[:branches],(all_newvar[JuMP.index(branch[1]).value],branch[2],branch[3]))
		elseif typeof(branch[1])==GenericAffExpr{Float64,VariableRef}
			exp=0
			for (var,coef) in branch[1].terms
				exp+=coef*all_newvar[JuMP.index(var).value]
			end
			push!(newmodel.master_problem.ext[:branches],(exp,branch[2],branch[3]))
		end
		add_branch_constraint(newmodel.master_problem,newmodel.master_problem.ext[:branches][end])
	end

	if typeof(branch[1])==VariableRef
		branch=(all_newvar[JuMP.index(branch[1]).value],branch[2],branch[3])
	elseif typeof(branch[1])==GenericAffExpr{Float64,VariableRef}
		exp=0
		for (var,coef) in branch[1].terms
			exp+=coef*all_newvar[JuMP.index(var).value]
		end
		branch=(exp,branch[2],branch[3])
	end
	new_branch_constraint(newmodel.master_problem,branch)

	nodes=collect(jmodel.tree)
	for i in num+1:length(all_var)
		for node in nodes
			if normalized_coefficient(jmodel.master_problem.ext[:convexcombination][node], all_var[i]) == 1.0
				column=copy_column(jmodel.master_problem, jmodel.sub_problems[node], node, all_var[i])
				newvar=add_variable_as_column(newmodel.master_problem, UnitIntervalInformation(), column)
				break
			end
		end
	end

	return newmodel
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

#function variable_branch(master::JuMP.Model, tree::T where T <: AbstractTree, expansions::Dict{AbstractTree,Dict{Symbol,Any}}, inttol::Float64)
function variable_branch(master, tree, expansions, inttol)
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

function constraint_branch(master, tree, expansions, inttol)
	branches=Array{Any,1}()
	parents=history(tree)
	biggest_frac=inttol
	for leaf in collect(tree)
		if typeof(leaf)==Leaf
			for x in keys(expansions[leaf])
				var = expansions[leaf][x]
				if typeof(var) <: AbstractArray
					for key in keys(var)
						nodes=[]
						for node in parents(leaf)
							if min(JuMP.value(expansions[node][x][key]),1-JuMP.value(expansions[node][x][key]))>inttol
								push!(nodes,node)
							end
						end
						if length(nodes)!=0
							for node in nodes
								expr=@expression(master,expansions[node][x][key])
								push!(branches,(expr,:eq,1))
							end
							expr=@expression(master,sum(expansions[n][x][key] for n in nodes))
							push!(branches,(expr,:eq,0))
							return branches
						end
					end
				else
					nodes=[]
					for node in parents(leaf)
						if min(JuMP.value(expansions[node][x]),1-JuMP.value(expansions[node][x]))>inttol
							push!(nodes,node)
						end
					end
					if length(nodes)!=0
						for node in nodes
							expr=@expression(master,expansions[node][x])
							push!(branches,(expr,:eq,1))
						end
						expr=@expression(master,sum(expansions[n][x] for n in nodes))
						push!(branches,(expr,:eq,0))
						return branches
					end
				end
			end
		end
	end
	branches
end

function perform_branch(jmodel::JuDGEModel,branches::Array{Any,1})
	newmodels=Array{JuDGEModel,1}()

	for branch in branches
		push!(newmodels,copy_model(jmodel,branch))
	end

	return newmodels
end

function branch_and_price(judge::JuDGEModel;branch_method=JuDGE.constraint_branch,search=:depth_first_resurface,
   abstol=10^-14,
   reltol=10^-14,
   rlx_abstol=10^-14,
   rlx_reltol=10^-14,
   duration=Inf,
   iter=2^63-1,
   inttol=10^-8,
   column_subset=1.0)

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

	judge.master_problem.ext[:branches]=Array{Any,1}()
	judge.master_problem.ext[:newbranches]=Array{Any,1}()

	models = Array{JuDGEModel,1}()
	newmodels = Array{JuDGEModel,1}()

	push!(models,judge)

	UB=Inf
	i=1
	best=0

	while true
		model=models[i]
		LB=Inf
		bestLB=0

		N=length(models)
		while model.bounds.LB>UB
		  if i==N
      	  	return models[best]
		  end
		  i+=1
		  model=models[i]
	   end

		for j in i:N
			if models[j].bounds.LB<LB
			 	LB=models[j].bounds.LB
				bestLB=j
			end
		end

		if UB-LB<abstol || (UB-LB)/(abs(UB))<reltol
			println("Integer solution meets convergence tolerances")
			return models[best]
		end

		if search==:lowestLB
			model=models[bestLB]
			deleteat!(models,bestLB)
			insert!(models,i,model)
		end
		println("Model "*string(i)*" of "*string(N)*". UB: "*string(UB)*", LB:"*string(LB))

		if rlx_reltol_func
			rrt=rlx_reltol(i,N)
		end

		if rlx_abstol_func
			rat=rlx_abstol(i,N)
		end

		solve(model,abstol=abstol,reltol=reltol,rlx_abstol=rat,rlx_reltol=rrt,duration=duration,iter=iter,inttol=inttol,allow_frac=1,prune=UB,column_subset=column_subset)

		if model.bounds.UB<UB
			UB=model.bounds.UB
			best=i
		end

		if model.bounds.LB<=UB && termination_status(model.master_problem)==MathOptInterface.OPTIMAL
			branches=branch_method(model.master_problem, model.tree, model.master_problem.ext[:expansions], inttol)
			if length(branches)>0
			    newmodels=perform_branch(model,branches)
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
				return models[best]
			else
				println("Found integer solution.")
				i+=1
			end
		elseif i==length(models)
			solve_binary(models[best])
			println("Objective value of best integer-feasible solution: "*string(objective_value(models[best].master_problem)))
			return models[best]
		else
			i+=1
		end
	end


end
