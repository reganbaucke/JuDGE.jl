module JuDGE

using JuMP
using MathOptInterface
using Printf

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

struct JuDGEModel
   tree::AbstractTree
   master_problem::JuMP.Model
   sub_problems::Dict{AbstractTree,JuMP.Model}
	bounds::Bounds
	discount_factor::Float64
	master_solver
	# function JuDGEModel(tree::T where T <: AbstractTree,master_problem::JuMP.Model,sub_problems::Dict{AbstractTree,JuMP.Model},bounds::Bounds,discount_factor::Float64)
	# 	return new(tree,master_problem,sub_problems,bounds,discount_factor)
	# end

	function JuDGEModel(tree::T where T <: AbstractTree, probability_function, sub_problems, solver, bounds, discount_factor::Float64)
 	  probabilities = probability_function(tree)
       master_problem = build_master(sub_problems, tree, probabilities, solver, discount_factor)
       return new(tree,master_problem,sub_problems,Bounds(bounds.UB,bounds.LB),discount_factor,solver)
    end

   function JuDGEModel(tree::T where T <: AbstractTree, probability_function, sub_problem_builder, solver; discount_factor=1.0)
      #this = new()
      println("Establishing JuDGE model for tree: " * string(tree))
	  probabilities = probability_function(tree)
	  sub_problems = Dict(i => sub_problem_builder(i) for i in collect(tree))
      scale_objectives(tree,sub_problems,probabilities,discount_factor)
      print("Checking sub-problem format...")
      check_specification_is_legal(sub_problems)
      println("Passed")
      print("Building master problem...")
      master_problem = build_master(sub_problems, tree, probabilities, solver, discount_factor)
      println("Complete")
      return new(tree,master_problem,sub_problems, Bounds(Inf,-Inf),discount_factor,solver)
   end
end

include("branchandprice.jl")

function build_master(sub_problems, tree::T where T <: AbstractTree, probabilities, solver, discount_factor::Float64)
   model = Model(solver)
   @objective(model,Min,0)

   history_function = history(tree)
   depth_function = depth(tree)

   # load in the variables
   model.ext[:expansions] = Dict{AbstractTree,Dict{Symbol,Any}}()
   for (node,sp) in sub_problems
      model.ext[:expansions][node] = Dict{Symbol,Any}()
      for (name,variable) in sp.ext[:expansions]
         model.ext[:expansions][node][name] = copy_variable!(model, variable, relaxbinary)
      end
   end

   # do the object function for the master
   # should be able to implement this with for each
   for (node,sp) in sub_problems
	  df=discount_factor^depth_function(node)
      for (name,variable) in sp.ext[:expansions]
         if typeof(variable) <: AbstractArray
            for i in eachindex(variable)
               set_objective_coefficient(model, model.ext[:expansions][node][name][i], df*probabilities(node)*coef(sp.ext[:expansioncosts],variable[i]))
            end
         else
            set_objective_coefficient(model, model.ext[:expansions][node][name], df*probabilities(node)*coef(sp.ext[:expansioncosts],variable))
         end
      end
   end

   # create the cover constraints
   model.ext[:coverconstraint] = Dict{AbstractTree,Dict{Symbol,Any}}()
   for (node,sp) in sub_problems
      model.ext[:coverconstraint][node] = Dict{Symbol,Any}()
      for (name,variable) in sp.ext[:expansions]
         if typeof(variable) <: AbstractArray
            model.ext[:coverconstraint][node][name] = Dict()
            for i in eachindex(variable)
               model.ext[:coverconstraint][node][name][i] = @constraint(model, 0 <= sum(model.ext[:expansions][past][name][i] for past in history_function(node)))
            end
			if typeof(node)==Leaf
				for i in eachindex(variable)
	            	@constraint(model, sum(model.ext[:expansions][past][name][i] for past in history_function(node))<=1)
	            end
			end
         else
            model.ext[:coverconstraint][node][name] = @constraint(model, 0 <= sum(model.ext[:expansions][past][name] for past in history_function(node)))
			if typeof(node)==Leaf
				@constraint(model, sum(model.ext[:expansions][past][name] for past in history_function(node))<=1)
			end
         end
      end
   end

   model.ext[:convexcombination] = Dict{AbstractTree,ConstraintRef}()
   for node in keys(sub_problems)
      model.ext[:convexcombination][node] = @constraint(model, 0 == 1)
   end

   model
end


function scale_objectives(tree::T where T <: AbstractTree,sub_problems, probabilities,discount_factor::Float64)
   depth_function=depth(tree)

   for (node,sp) in sub_problems
      @objective(sp, Min, objective_function(sp)*probabilities(node)*(discount_factor^depth_function(node)))
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

	set_objective_coefficient(master, holder, column.obj)
end

function build_column(master, sub_problem ,node)
   singlevars=Array{Any,1}()
   arrayvars=Dict{Symbol,Array{Any,1}}()

   for (name,variable) in sub_problem.ext[:expansions]
      if typeof(variable) <: AbstractArray
		 arrayvars[name]=Array{Int64,1}()
         for i in eachindex(variable)
            if JuMP.value(variable[i]) >= 0.99
               push!(arrayvars[name],i)
            end
         end
      else
         if JuMP.value(variable) >= 0.99
            push!(singlevars,name)
         end
      end
   end
   Column(node,singlevars,arrayvars,get_objective_coef_for_column(sub_problem))
end

function get_objective_coef_for_column(sub_problem)
   ### objective coefficient is the objective function value minus the terms with expansions
   coef = objective_value(sub_problem)
   for (name,var) in sub_problem.ext[:expansions]
      if typeof(var) <: AbstractArray
         for single_var in var
            coef -= JuMP.value(single_var)*objcoef(single_var)
         end
      else
         coef -= JuMP.value(var)*objcoef(var)
      end
   end
   return coef
end

function solve(judge::JuDGEModel;
   abstol= 10^-14,
   reltol= 10^-14,
   rlx_abstol= 10^-14,
   rlx_reltol= 10^-14,
   duration= Inf,
   iter= 2^63 - 1,
   inttol=10^-9, # The Maximum int
   allow_frac=0,
   prune=Inf,
   )

   # encode the user convergence test in a ConvergenceState struct
   done = ConvergenceState(0.0, 0.0, 0.0, abstol, reltol, rlx_abstol, rlx_reltol, duration, iter, inttol)

   current = InitialConvergenceState()

	print("")
	print("")
  	println("Current ObjVal  |   Upper Bound   Lower Bound  |  Absolute Diff   Relative Diff  |  Fractionality  |      Time     Iter")
   # set up times for use in convergence
   initial_time = time()
   stamp = initial_time
   obj = Inf
   while true
      # perform the main iterations
      optimize!(judge.master_problem)
      status=termination_status(judge.master_problem)

      for node in collect(judge.tree)
         updateduals(judge.master_problem, judge.sub_problems[node],node, status)
         optimize!(judge.sub_problems[node])
      end

	  frac=NaN
      if status==MathOptInterface.OPTIMAL
		  getlowerbound(judge)
		  frac = absolutefractionality(judge)
	      obj = objective_value(judge.master_problem)
	  	  if frac<done.int
			  judge.bounds.UB=obj
		  end
	  elseif judge.bounds.LB>-Inf
		  println("\nMaster problem is infeasible")
		  return
	  end

      current = ConvergenceState(obj, judge.bounds.UB, judge.bounds.LB, time() - initial_time, current.iter + 1, frac)
      println(current)

	  if has_converged(done, current) || prune<judge.bounds.LB || (allow_frac==2 && frac>done.int)
		  break
	  end

	  for node in collect(judge.tree)
	  	  column = build_column(judge.master_problem, judge.sub_problems[node], node)
  		  add_variable_as_column(judge.master_problem, UnitIntervalInformation(), column)
	  end
   end

   if allow_frac==0 && current.int>done.int
		solve_binary(judge)
		current = ConvergenceState(judge.bounds.UB, judge.bounds.UB, judge.bounds.LB, time() - initial_time, current.iter + 1, absolutefractionality(judge))
		println(current)
   end
   println("\nConvergence criteria met.")
end

function getlowerbound(judge::JuDGEModel)
   lb = objective_value(judge.master_problem)
   for (node,sp) in judge.sub_problems
      lb += objective_value(sp)-dual(judge.master_problem.ext[:convexcombination][node])
   end
   if lb>judge.bounds.LB
	  judge.bounds.LB=lb
   end
end

function solve_binary(judge::JuDGEModel)
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
   optimize!(judge.master_problem)
   judge.bounds.UB=objective_value(judge.master_problem)

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

function updateduals(master, sub_problem, node, status)
   for (name,var) in sub_problem.ext[:expansions]
      if typeof(var) <: AbstractArray
         for i in eachindex(var)
            if status == MathOptInterface.OPTIMAL
               set_objective_coefficient(sub_problem, var[i], -dual(master.ext[:coverconstraint][node][name][i]))
            else
               set_objective_coefficient(sub_problem, var[i], -9999.0)
            end
         end
      else
         if status == MathOptInterface.OPTIMAL
            set_objective_coefficient(sub_problem, var, -dual(master.ext[:coverconstraint][node][name]))
         else
            set_objective_coefficient(sub_problem, var, -9999.0)
         end
      end
   end
end

function fix_expansions(jmodel::JuDGEModel;node=jmodel.tree::AbstractTree,invest0=Dict{Any,Float64}())
   if termination_status(jmodel.master_problem) != MathOptInterface.OPTIMAL
      error("You need to first solve the decomposed model.")
   end

   # if jmodel.is_fixed && node==jmodel.tree
   #    error("You have already fixed this model previously, you need to rebuid it.")
   # end

   invest=copy(invest0)

   queue=[]
   for key in keys(jmodel.master_problem.ext[:expansions][node])
      var = jmodel.master_problem.ext[:expansions][node][key]
      var2 = jmodel.sub_problems[node][key]
      if isa(var,Array)
         for v in keys(var)
            push!(queue,v)
            if v in keys(invest)
               invest[v]+=JuMP.value(var[v])
            else
               invest[v]=JuMP.value(var[v])
            end
         end
         index=1

         for v in var2
            LHS=AffExpr(0.0)
            add_to_expression!(LHS,1.0,v)
			if invest[queue[index]]>0.99
				@constraint(jmodel.sub_problems[node],LHS<=1.0)
			else
				@constraint(jmodel.sub_problems[node],LHS<=0.0)
			end
            @constraint(jmodel.sub_problems[node],LHS<=invest[queue[index]])
            set_objective_coefficient(jmodel.sub_problems[node], v, 0.0)
            index+=1
         end
      elseif isa(var,VariableRef)
         if string(var2) in keys(invest)
            invest[string(var2)]+=JuMP.value(var)
         else
            invest[string(var2)]=JuMP.value(var)
         end

         LHS=AffExpr(0.0)
         add_to_expression!(LHS,1.0,var2)
         @constraint(jmodel.sub_problems[node],LHS<=invest[string(var2)])
         set_objective_coefficient(jmodel.sub_problems[node], var2, 0.0)
      elseif isa(var,JuMP.Containers.DenseAxisArray) || isa(var,JuMP.Containers.SparseAxisArray)
         val=JuMP.value.(var)
         for key2 in keys(var)
            if (key,key2) in keys(invest)
               invest[(key,key2)]+=val[key2]
            else
               invest[(key,key2)]=val[key2]
            end
            LHS=AffExpr(0.0)
            add_to_expression!(LHS,1.0,var2[key2])
            @constraint(jmodel.sub_problems[node],LHS<=invest[(key,key2)])
            set_objective_coefficient(jmodel.sub_problems[node], var2[key2], 0.0)
        end
      end
   end

   if typeof(node)==Tree
      for i in 1:length(node.children)
         fix_expansions(jmodel,node=node.children[i],invest0=invest)
      end
   end
end

function resolve_fixed(jmodel::JuDGEModel)
   obj=0.0
   for n in collect(jmodel.tree)
      JuMP.optimize!(jmodel.sub_problems[n])
      obj+=objective_value(jmodel.sub_problems[n])
      for key in keys(jmodel.master_problem.ext[:expansions][n])
         var = jmodel.master_problem.ext[:expansions][n][key]
         if isa(var,Array)
            for v in keys(var)
               obj+=JuMP.value(var[v])*objcoef(var[v])
            end
         elseif isa(var,VariableRef)
            obj+=JuMP.value(var)*objcoef(var)
         elseif isa(var,JuMP.Containers.DenseAxisArray) || isa(var,JuMP.Containers.SparseAxisArray)
            val=JuMP.value.(var)
            for key2 in keys(val)
               obj+=val[key2]*objcoef(var[key2])
            end
         end
      end
    end
    return obj
end

include("model_verification.jl")

# pretty printing
function Base.show(io::IO, ::MIME"text/plain", judge::JuDGEModel)
   print(io, "JuDGE Model with:\n")
   println(io, "  Tree: ", judge.tree)
   print(io, "  Expansion variables: ")
   keys = collect(get_expansion_keys(judge.sub_problems[judge.tree]))
   for i in 1:length(keys)
      print(io, "$(keys[i]) ")
   end
end

function Base.show(io::IO, judge::JuDGEModel)
   print(io, "JuDGE Model with:\n")
   println(io, "  Tree: ", judge.tree)
   print(io, "  Expansion variables: ")
   keys = collect(get_expansion_keys(judge.sub_problems[judge.tree]))
   for i in 1:length(keys)
      print(io, "$(keys[i]) ")
   end
end
include("output.jl")

export @expansion, @expansionconstraint, @expansioncosts, JuDGEModel, Leaf, Tree, AbstractTree, narytree, ConditionallyUniformProbabilities, get_node, tree_from_leaves, tree_from_nodes, tree_from_file, DetEqModel

end
