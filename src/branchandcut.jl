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

function new_branch_constraint(master::JuMP.Model,index::Int64)
	add_branch_constraint(master,master.ext[:newbranches][index])

	push!(master.ext[:branches],master.ext[:newbranches][index])
	while length(master.ext[:newbranches])>0
		deleteat!(master.ext[:newbranches],1)
	end
end

function copy_model(jmodel::JuDGEModel)
	newmodel=JuDGEModel(jmodel.tree, ConditionallyUniformProbabilities, jmodel.sub_problems, jmodel.master_solver, Bounds(Inf,jmodel.bounds.LB), jmodel.discount_factor)

	newmodel.master_problem.ext[:branches]=Array{Any,1}()
	newmodel.master_problem.ext[:newbranches]=Array{Any,1}()

	all_newvar=all_variables(newmodel.master_problem)
	all_var=all_variables(jmodel.master_problem)
	num=length(all_newvar)

	for branch in jmodel.master_problem.ext[:branches]
		if typeof(branch[1])==VariableRef
			push!(newmodel.master_problem.ext[:branches],(all_newvar[JuMP.index(branch[1]).value],branch[2],branch[3]))
		elseif typeof(branch[1])==GenericAffExpr{Float64,VariableRef}
			exp=0
			for (coef,var) in branch[1].terms
				exp+=coef*all_newvar[JuMP.index(var).value]
			end
			push!(newmodel.master_problem.ext[:branches],(exp,branch[2],branch[3]))
		end
		add_branch_constraint(newmodel.master_problem,newmodel.master_problem.ext[:branches][end])
	end

	for branch in jmodel.master_problem.ext[:newbranches]
		if typeof(branch[1])==VariableRef
			push!(newmodel.master_problem.ext[:newbranches],(all_newvar[JuMP.index(branch[1]).value],branch[2],branch[3]))
		elseif typeof(branch[1])==GenericAffExpr{Float64,VariableRef}
			exp=0
			for (coef,var) in branch[1].terms
				exp+=coef*all_newvar[JuMP.index(var).value]
			end
			push!(newmodel.master_problem.ext[:newbranches],(exp,branch[2],branch[3]))
		end
	end

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

function most_fractional(jmodel::JuDGEModel)
	branch=nothing
	maxfrac=0.000001
	for node in collect(jmodel.tree)
		for x in keys(jmodel.master_problem.ext[:expansions][node])
			var = jmodel.master_problem.ext[:expansions][node][x]
			if typeof(var) <: AbstractArray
				for key in keys(var)
					fractionality=min(JuMP.value(var[key]),1-JuMP.value(var[key]))
					if fractionality>maxfrac
						branch=@expression(jmodel.master_problem,var[key])
						maxfrac=fractionality
					end
				end
			else
			  fractionality=min(JuMP.value(var),1-JuMP.value(var))
			  if fractionality>maxfrac
				  branch=@expression(jmodel.master_problem,var)
				  maxfrac=fractionality
			  end
			end
		end
	end
	if branch==nothing
		return false
	end
	push!(jmodel.master_problem.ext[:newbranches],(branch,:eq,0))
	push!(jmodel.master_problem.ext[:newbranches],(branch,:eq,1))
	return true
end

function first_fraction(jmodel::JuDGEModel;inttol=10^-9)
	branch=nothing
	parents=history(jmodel.tree)
	biggest_frac=inttol
	for node in collect(jmodel.tree)
		if typeof(node)==Leaf
			continue
		end
		for x in keys(jmodel.master_problem.ext[:expansions][node])
			var = jmodel.master_problem.ext[:expansions][node][x]
			if typeof(var) <: AbstractArray
				for key in keys(var)
					if min(JuMP.value(var[key]),1-JuMP.value(var[key]))>biggest_frac
						biggest_frac=min(JuMP.value(var[key]),1-JuMP.value(var[key]))
						branch=@expression(jmodel.master_problem,sum(jmodel.master_problem.ext[:expansions][n][x][key] for n in parents(node)))
						push!(jmodel.master_problem.ext[:newbranches],(branch,:eq,1))
						push!(jmodel.master_problem.ext[:newbranches],(branch,:eq,0))
						#
						# for n in collect(node)
						# 	if n!=node
						# 		push!(vars1,jmodel.master_problem.ext[:expansions][n][x][key])
						# 	end
						# end
						return true
					end
				end
			elseif min(JuMP.value(var),1-JuMP.value(var))>biggest_frac
				biggest_frac=min(JuMP.value(var[key]),1-JuMP.value(var))
				branch=@expression(jmodel.master_problem,sum(jmodel.master_problem.ext[:expansions][n][x] for n in parents(node)))
				push!(jmodel.master_problem.ext[:newbranches],(branch,:eq,0))
				push!(jmodel.master_problem.ext[:newbranches],(branch,:eq,1))
					# for n in collect(node)
					# 	if n!=node
					# 		push!(vars1,jmodel.master_problem.ext[:expansions][n][x][key])
					# 	end
					# end
				  return true
			end
		end
	end
	return false
end

function perform_branch(jmodel::JuDGEModel)
	newmodels=Array{JuDGEModel,1}()

	push!(newmodels,jmodel)
	for i in 2:(length(jmodel.master_problem.ext[:newbranches]))
		push!(newmodels,copy_model(jmodel))
	end

	for i in 1:length(newmodels)
		new_branch_constraint(newmodels[i].master_problem,i)
	end

	return newmodels
end

function branch_and_price(judge::JuDGEModel;branch_method=JuDGE.most_fractional,search=:depth_first_resurface,
   abstol= 10^-14,
   reltol= 10^-14,
   rlx_abstol= 10^-14,
   rlx_reltol= 10^-14,
   duration= Inf,
   iter= 2^63 - 1,
   inttol=10^-14)

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

		while model.bounds.LB>UB
		  if i==length(models)
      	  	return models[best]
		  end
		  i=i+1
		  model=models[i]
	   end

		for j in i:length(models)
			if models[j].bounds.LB<LB
			 	LB=models[j].bounds.LB
				bestLB=j
			end
		end

		if search==:lowestLB
			temp=models[bestLB]
			deleteat!(models,bestLB)
			insert!(models,i,temp)
		end
		println("Model "*string(i)*" of "*string(length(models))*". UB: "*string(UB)*", LB:"*string(LB))
		solve(model,abstol=abstol,reltol=reltol,rlx_abstol=max(rlx_abstol,0.01/length(models)),rlx_reltol=rlx_reltol,duration=duration,iter=iter,inttol=inttol,allow_frac=1,prune=UB)
		if model.bounds.UB<UB
			UB=model.bounds.UB
			best=i
		end

		if model.bounds.LB<=UB
			if most_fractional(model)
			#if first_fraction(model)
			   newmodels=perform_branch(model)
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
				i+=1
			end
		elseif i==length(models)
			return models[best]
		else
			i+=1
		end
	end
end
