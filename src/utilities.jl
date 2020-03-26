# helper function which fetches the current objective value coef for a given variable
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
