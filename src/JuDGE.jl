module JuDGE

using JuMP
using MathOptInterface
using Printf

include("tree.jl")
include("macros.jl")
include("convergence.jl")
include("utilities.jl")
include("deteq.jl")

struct JuDGEModel
   tree::AbstractTree
   master_problem::JuMP.Model
   sub_problems::Dict{AbstractTree,JuMP.Model}
   function JuDGEModel(tree, probability_function, sub_problem_builder, solver)
      #this = new()
      println("Establishing JuDGE model for tree: " * string(tree))
      probabilities = probability_function(tree)
      sub_problems = Dict(i => sub_problem_builder(i) for i in collect(tree))
      scale_objectives(sub_problems,probabilities)
      print("Checking sub-problem format...")
      check_specification_is_legal(sub_problems)
      println("Passed")
      print("Building master problem...")
      master_problem = build_master(sub_problems, tree, probabilities, solver)
      println("Complete")
      return new(tree,master_problem,sub_problems)
   end
end

function build_master(sub_problems, tree::T where T <: AbstractTree, probabilities, solver)
   model = Model(solver)
   @objective(model,Min,0)

   history_function = history(tree)

   # load in the variables
   model.ext[:expansions] = Dict()
   for node in keys(sub_problems)
      model.ext[:expansions][node] = Dict()
      for variable in sub_problems[node].ext[:expansions]
         name = get_variable_name(sub_problems[node], variable)
         model.ext[:expansions][node][name] = copy_variable!(model, variable, relaxbinary)
      end
   end

   # do the object function for the master
   # should be able to implement this with for each
   for node in keys(sub_problems)
      for variable in sub_problems[node].ext[:expansions]
         if typeof(variable) <: AbstractArray
            name = get_variable_name(sub_problems[node], variable)
            for i in eachindex(variable)
               set_objective_coefficient(model, model.ext[:expansions][node][name][i], probabilities(node)*coef(sub_problems[node].ext[:expansioncosts],variable[i]))
            end
         else
            name = get_variable_name(sub_problems[node], variable)
            set_objective_coefficient(model, model.ext[:expansions][node][name], probabilities(node)*coef(sub_problems[node].ext[:expansioncosts],variable))
         end
      end
   end

   # create the cover constraints
   model.ext[:coverconstraint] = Dict()
   for node in keys(sub_problems)
      model.ext[:coverconstraint][node] = Dict()
      for variable in sub_problems[node].ext[:expansions]
         if typeof(variable) <: AbstractArray
            name = get_variable_name(sub_problems[node], variable)
            model.ext[:coverconstraint][node][name] = Dict()
            for i in eachindex(variable)
               model.ext[:coverconstraint][node][name][i] = @constraint(model, 0 <= sum(model.ext[:expansions][past][name][i] for past in history_function(node)))
            end
         else
            name = get_variable_name(sub_problems[node], variable)
            model.ext[:coverconstraint][node][name] = @constraint(model, 0 <= sum(model.ext[:expansions][past][name] for past in history_function(node)))
         end
      end
   end

   model.ext[:convexcombination] = Dict()
   for node in keys(sub_problems)
      model.ext[:convexcombination][node] = @constraint(model, 0 == 1)
   end

   model
end


function scale_objectives(sub_problems, probabilities)
   for node in keys(sub_problems)
      @objective(sub_problems[node], Min, objective_function(sub_problems[node])*probabilities(node))
   end
end

function Base.map(f, hello::JuMP.Containers.DenseAxisArray)
   JuMP.Containers.DenseAxisArray(map(f, hello.data), deepcopy(hello.axes),deepcopy(hello.lookup))
end

function Base.map(f, hello::JuMP.Containers.SparseAxisArray)
   JuMP.Containers.SparseAxisArray(Dict(i => f(hello[i]) for i in keys(hello.data)))
end

function add_variable_as_column(model, info, objcoef, constraints)
   holder = JuMP.add_variable(model, JuMP.build_variable(error, info))
   foreach(constraints) do con
      set_normalized_coefficient(con, holder, 1.0)
   end
   set_objective_coefficient(model, holder, objcoef)
end

function build_column(master, sub_problem, node)
   (get_objective_coef_for_column(sub_problem), collect_constraints(master,sub_problem, node))
end

function collect_constraints(master, sub_problem ,node)
   collection = []
   push!(collection, master.ext[:convexcombination][node])
   for variable in sub_problem.ext[:expansions]
      if typeof(variable) <: AbstractArray
         for i in eachindex(variable)
            name = get_variable_name(sub_problem, variable)
            if JuMP.value(variable[i]) == 1.0
               push!(collection, master.ext[:coverconstraint][node][name][i])
            end
         end
      else
         if JuMP.value(variable) == 1.0
            name = get_variable_name(sub_problem, variable)
            push!(collection, master.ext[:coverconstraint][node][name])
         end
      end
   end
   collection
end

function get_objective_coef_for_column(sub_problem)
   ### objective coefficient is the objective function value minus the terms with expansions
   coef = objective_value(sub_problem)
   for var in sub_problem.ext[:expansions]
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
   abstol= 0,
   reltol= -Inf,
   duration= Inf,
   iter= 2^63 - 1) # The Maximum int

   # encode the user convergence test in a ConvergenceState struct
   done = ConvergenceState(0, 0, abstol, reltol, duration, iter)

   current = InitialConvergenceState()

   Printf.@printf("\n    Upper Bound     Lower Bound   Absolute Diff   Relative Diff       Time    Iter\n")
   println(current)

   # set up times for use in convergence
   initial_time = time()
   stamp = initial_time

   LB=-Inf

   while !has_converged(done, current)
      # perform the main iterations
      JuMP.optimize!(judge.master_problem)
      status=termination_status(judge.master_problem)
      upper_bound = objective_value(judge.master_problem)
      for node in collect(judge.tree)
         updateduals(judge.master_problem, judge.sub_problems[node],node, status)
         optimize!(judge.sub_problems[node])
         (obj_coef, constraints) = build_column(judge.master_problem, judge.sub_problems[node], node)
         add_variable_as_column(judge.master_problem, UnitIntervalInformation(), obj_coef, constraints)
      end

      # update current convergence state
      lb = getlowerbound(judge)
      if lb>LB
         LB=lb
      else
         lb=LB
      end
      ub = objective_value(judge.master_problem)
      current = ConvergenceState(ub, lb, ub - lb, (ub - lb)/abs(ub), time() - initial_time, current.iter + 1)
      #### Print a small update on the convergence stats
      if time() - stamp > 0.0
         stamp = time()
         println(current)
      end
   end
   JuMP.optimize!(judge.master_problem)

   println("\nConvergence criteria met.")

   judge
end

function getlowerbound(judge::JuDGEModel)
   JuMP.optimize!(judge.master_problem)
   lb = objective_value(judge.master_problem)
   for n in keys(judge.sub_problems)
      lb += objective_value(judge.sub_problems[n]) - dual(judge.master_problem.ext[:convexcombination][n])
   end
   return lb
end

function updateduals(master, sub_problem, node, status)
   for var in sub_problem.ext[:expansions]
      if typeof(var) <: AbstractArray
         name = get_variable_name(sub_problem,var)
         for i in eachindex(var)
            if status == MathOptInterface.OPTIMAL
               optimize!(master)
               set_objective_coefficient(sub_problem, var[i], -dual(master.ext[:coverconstraint][node][name][i]))
            else
               # set_objective_coefficient(sub_problem, var[i], -dual(master.ext[:coverconstraint][node][name][i]))
               set_objective_coefficient(sub_problem, var[i], -9999.0)
            end
         end
      else
         name = get_variable_name(sub_problem,var)
         if status == MathOptInterface.OPTIMAL
            optimize!(master)
            set_objective_coefficient(sub_problem, var, -dual(master.ext[:coverconstraint][node][name]))
            #set_objective_coefficient(sub_problem, var, -dual(master.ext[:coverconstraint][node][name]))
         else
            # set_objective_coefficient(sub_problem, var, -dual(master.ext[:coverconstraint][node][name]))
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
               invest[v]+=JuMP.value(var[v...])
            else
               invest[v]=JuMP.value(var[v...])
            end
         end
         index=1

         for v in var2
            LHS=AffExpr(0.0)
            add_to_expression!(LHS,1.0,v)
            @constraint(jmodel.sub_problems[node],LHS==invest[queue[index]])
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
         @constraint(jmodel.sub_problems[node],LHS==invest[string(var2)])
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
            @constraint(jmodel.sub_problems[node],LHS==invest[(key,key2)])
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
               obj+=JuMP.value(var[v...])*objcoef(var[v...])
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
   for i in 1:length(keys)-1
      print(io, "$(keys[i]), ")
   end
   print(io, "$(keys[end]).")
end

function Base.show(io::IO, judge::JuDGEModel)
   print(io, "JuDGE Model with:\n")
   println(io, "  Tree: ", judge.tree)
   print(io, "  Expansion variables: ")
   keys = collect(get_expansion_keys(judge.sub_problems[judge.tree]))
   for i in 1:length(keys)-1
      print(io, "$(keys[i]), ")
   end
   print(io, "$(keys[end]).")
end
include("output.jl")

export @expansion, @expansionconstraint, @expansioncosts, JuDGEModel, Leaf, Tree, AbstractTree, narytree, ConditionallyUniformProbabilities, get_node, tree_from_leaves, tree_from_nodes, tree_from_file, DetEqModel

end
