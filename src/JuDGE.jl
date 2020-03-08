module JuDGE

using JuMP
using MathOptInterface

include("tree.jl")
include("macros.jl")
include("convergence.jl")


function copy_variable!(toModel, variable)
   map(variable) do x
      JuMP.add_variable(toModel, JuMP.build_variable(error, (get_info(x))))
   end
end

function copy_variable!(toModel, variable, f)
   map(variable) do x
      JuMP.add_variable(toModel, JuMP.build_variable(error, f(get_info(x))))
   end
end

function copy_variable!(toModel, variable::JuMP.VariableRef)
   JuMP.add_variable(toModel, JuMP.build_variable(error, get_info(variable)))
end

function copy_variable!(toModel, variable::JuMP.VariableRef, f)
   JuMP.add_variable(toModel, JuMP.build_variable(error, f(get_info(variable))), "exp")
end


# constuct variable info object for a single variable
function get_info(x::VariableRef)
   has_lb_local = false
   lb_local = NaN
   has_ub_local = false
   ub_local = NaN
   is_fixed_local = false
   fixed_value_local = NaN
   has_start_local = false
   start_value_local = NaN
   is_binary_local = false
   is_integer_local = false

   if has_lower_bound(x)
      has_lb_local = true
      lb_local = lower_bound(x)
   end

   if has_upper_bound(x)
      has_ub_local = true
      ub_local = upper_bound(x)
   end

   if is_fixed(x)
      is_fixed_local = true
      fixed_value_local = fix_value(x)
   end

   if start_value(x) != nothing
      has_start_local = true
      start_value_local = start_value(x)
   end

   if is_binary(x)
      is_binary_local = true
   end

   if is_integer(x)
      is_integer_local = true
   end

   VariableInfo(has_lb_local, lb_local, has_ub_local, ub_local, is_fixed_local, fixed_value_local, has_start_local, start_value_local, is_binary_local, is_integer_local)
end

function relaxbinary(x::VariableInfo)
   VariableInfo(true, 0.0, true, 1.0, x.has_fix, x.fixed_value, x.has_start, x.start, false, x.integer)
end

function UnitIntervalInformation()
   VariableInfo(true, 0.0 , true, 1.0, false, NaN, false, NaN, false, false)
end

function objcoef(x::JuMP.VariableRef)
   affine_expression = objective_function(owner_model(x))
   if x in keys(affine_expression.terms)
      affine_expression.terms[x]
   else
      0.0
   end
end

function coef(aff, x::JuMP.VariableRef)
   if x in keys(aff.terms)
      aff.terms[x]
   else
      0.0
   end
end

function get_variable_name(sub_problem, variable)
   set = filter(keys(sub_problem.obj_dict)) do key
      if sub_problem.obj_dict[key] == variable
         true
      else
         false
      end
   end
   collect(set)[1]
end

struct JuDGEModel
   tree::AbstractTree
   master_problem::JuMP.Model
   sub_problems::Dict{AbstractTree,JuMP.Model}
   function JuDGEModel(tree, probability_function, sub_problem_builder, solver)
      #this = new()
      probabilities = probability_function(tree)
      sub_problems = Dict(i => sub_problem_builder(i) for i in collect(tree))
      scale_objectives(sub_problems,probabilities)
      check_specification_is_legal(sub_problems)
      master_problem = build_master(sub_problems, tree, probabilities, solver)
      return new(tree,master_problem,sub_problems)
   end
end

function build_master(sub_problems, tree::T where T <: AbstractTree, probabilities,solver)
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
            if value(variable[i]) == 1.0
               push!(collection, master.ext[:coverconstraint][node][name][i])
            end
         end
      else
         if value(variable) == 1.0
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
            coef -= value(single_var)*objcoef(single_var)
         end
      else
         coef -= value(var)*objcoef(var)
      end
   end
   return coef
end

function judgesolve(judge::JuDGEModel;
   abstol= -Inf,
   reltol= -Inf,
   duration= Inf,
   iter= 2^63 - 1) # The Maximum int

   # encode the user convergence test in a ConvergenceState struct
   done = ConvergenceState(abstol, reltol, duration, iter)
   # if the user put no options, then default to a convergence check on abstol
   if done == ConvergenceState(-Inf, -Inf, Inf, 2^63 - 1)
      done = ConvergenceState(0, -Inf, Inf, 2^63 - 1)
   end

   current = InitialConvergenceState()

   # set up times for use in convergence
   initial_time = time()
   stamp = initial_time

   LB=-Inf

   while !has_converged(done, current)
      # perform the main iterations
      JuMP.optimize!(judge.master_problem)

      upper_bound = objective_value(judge.master_problem)
      for node in collect(judge.tree)
         updateduals(judge.master_problem, judge.sub_problems[node],node, termination_status(judge.master_problem))
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
      current = ConvergenceState(ub - lb, (ub - lb)/abs(ub), time() - initial_time, current.iter + 1)
      #### Print a small update on the convergence stats
      if time() - stamp > 0.0
         stamp = time()
         println(current)
      end
   end
   JuMP.optimize!(judge.master_problem)

   println(current)

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

include("model_verification.jl")

# pretty printing
function Base.show(io::IO, ::MIME"text/plain", judge::JuDGEModel)
   print(io, "JuDGE Model with:\n")
   println(io, "  Tree: ", judge.tree)
   print(io, "  Expansion variables: ")
   keys = collect(get_expansion_keys(judge.sub_problems[judge.tree]))
   for i in 1:length(keys)-1
      print(io, "$(key[i]), ")
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

export @expansion, @expansionconstraint, @expansioncosts, JuDGEModel, judgesolve, history, Leaf, Tree, AbstractTree, narytree, ConditionallyUniformProbabilities, show, get_node, tree_from_leaves, tree_from_nodes, print_tree, JuDGE_value, print_expansions, tree_from_file

end
