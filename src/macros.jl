macro expansion(model, variable)
   ex = quote
      if !haskey(model.ext, :expansions)
         model.ext[:expansions] = Set()
      end
      push!(model.ext[:expansions], @variable(model, $variable))
   end
   return esc(ex)
end

macro expansionconstraint(model, variable)
   ex = quote
      if !haskey(model.ext, :expansionconstraints)
         model.ext[:expansionconstraints] = Set()
      end
      push!(model.ext[:expansionconstraints], @constraint(model, $variable))
   end
   return esc(ex)
end

macro expansioncosts(model, expr)
   ex = quote
      model.ext[:expansionconstraints] = @expression(model, $expr)
   end
   return esc(ex)
end

function copy_variable!(toModel, variable)
   new_copy = deepcopy(variable)

   for i in eachindex(variable)
      new_copy[i] = JuMP.add_variable(toModel, JuMP.build_variable(error, get_info(variable[i])))
   end

   new_copy
end

function copy_variable!(toModel, variable::JuMP.VariableRef)
   JuMP.add_variable(toModel, JuMP.build_variable(error, get_info(variable)))
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

optimize!(model2,with_optimizer(GLPK.Optimizer))
