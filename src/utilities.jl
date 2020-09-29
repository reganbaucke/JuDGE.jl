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

    VariableInfo(
        has_lb_local,
        lb_local,
        has_ub_local,
        ub_local,
        is_fixed_local,
        fixed_value_local,
        has_start_local,
        start_value_local,
        is_binary_local,
        is_integer_local,
    )
end

function relaxbinary(x::VariableInfo)
    VariableInfo(
        true,
        0.0,
        true,
        1.0,
        x.has_fix,
        x.fixed_value,
        x.has_start,
        x.start,
        false,
        x.integer,
    )
end

function UnitIntervalInformation()
    VariableInfo(true, 0.0, true, 1.0, false, NaN, false, NaN, false, false)
end

function objcoef(x::JuMP.VariableRef)
    affine_expression = objective_function(owner_model(x))
    if x in keys(affine_expression.terms)
        affine_expression.terms[x]
    else
        0.0
    end
end

# function coef(aff, x::JuMP.VariableRef)
#     if x in keys(aff.terms)
#         aff.terms[x]
#     else
#         0.0
#     end
# end

function coef(aff, x::JuMP.VariableRef)
    if x in keys(aff)
        aff[x]
    else
        0.0
    end
end

function unpack_expansions(a::Dict{AbstractTree,Dict{Symbol,Any}})
   assign=Array{Expr,1}()
   found=Dict{Symbol,Bool}()
   for (i,j) in a
       for (k,l) in j
           if k in keys(found)
               push!(assign,:($k[$i] = $l))
           else
               if isdefined(Main,k)
                   m=Symbol("→"*string(k))
                   push!(assign,:($m=$k))
               end
               push!(assign,:($k = Dict{AbstractTree,Any}()))
               push!(assign,:($k[$i] = $l))
               found[k]=true
           end
       end
   end
   assign
end

function clear_expansions(a::Dict{AbstractTree,Dict{Symbol,Any}})
   assign=Array{Expr,1}()
   for (i,j) in a
       for (k,l) in j
           m=Symbol("→"*string(k))
           if isdefined(Main,m)
               push!(assign,:($k=$m))
               push!(assign,:($m=nothing))
           else
               push!(assign,:($k=nothing))
           end
       end
       break
   end
   assign
end

function overprint(str)
    print("\u1b[0F")
    print("\u1b[0K")
    print(str)
end

function printleft(str)
    print("\u1b[40G")
    print("\u1b[1K")
    print("\u1b[1G")
    print(str)
end

function printright(str)
    print("\u1b[41G")
    print("\u1b[0K")
    print(str)
end
