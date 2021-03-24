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
        false,
    )
end

function relaxinteger(x::VariableInfo)
    VariableInfo(
        true,
        x.lower_bound,
        true,
        x.upper_bound,
        x.has_fix,
        x.fixed_value,
        x.has_start,
        x.start,
        false,
        false,
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

function densekey_to_tuple(key::Any)
    if typeof(key) <: JuMP.Containers.DenseAxisArrayKey
        return length(key.I)==1 ? key.I[1] : key.I
    else
        return key
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

function check_risk_type(risk::Any)
    passed=true
    if typeof(risk) <: Tuple{T,U} where T <: Union{Int64,Float64} where U <: Union{Int64, Float64}
        return
    elseif typeof(risk) <: Tuple{T,U,Dict{V,Float64}} where T <: Union{Int64,Float64} where U <: Union{Int64, Float64} where V <: AbstractTree
        return
    elseif typeof(risk) <: Array
        count=length(risk)
        for i in 1:length(risk)
            if typeof(risk[i]) <: Tuple{T,U} where T <: Union{Int64,Float64} where U <: Union{Int64, Float64}
                count-=1
            elseif typeof(risk[i]) <: Tuple{T,U,Dict{V,Float64}} where T <: Union{Int64,Float64} where U <: Union{Int64, Float64} where V <: AbstractTree
                count-=1
            else
                break
            end
        end
        if count==0
            return
        end
    end
    error("Invalid risk specification")
end

function compute_objval(scenarios::Dict{Leaf,Float64}, probabilities::Dict{AbstractTree,Float64}, risk)
	check_risk_type(risk)

    scenario_objs=Array{Tuple{Float64,Float64,Leaf},1}()
	EV_weight=1.0
	EV=0.0
	for (leaf,val) in scenarios
		push!(scenario_objs,(probabilities[leaf],val,leaf))
		EV+=probabilities[leaf]*val
	end

	obj=0.0

	if typeof(risk) <: Union{Tuple{Float64,Float64},Tuple{Float64,Float64,Dict{Leaf,Float64}}}
	    if risk[1]>0.0 && risk[2]<1.0
			EV_weight-=risk[1]
			so=scenario_objs
			if typeof(risk) <: Tuple{Float64,Float64,Dict{Leaf,Float64}}
				for j in 1:length(so)
					so[j]=(so[j][1],so[j][2]-risk[3][so[j][3]],so[j][3])
				end
			end
			sort!(so,by=i->i[2],rev=true)
			beta=risk[2]
			for scen in so
			 	if scen[1]>beta
			 		obj+=scen[2]*risk[1]*beta/risk[2]
			 		beta=0
			 	else
			 		obj+=scen[2]*risk[1]*scen[1]/risk[2]
			 		beta-=scen[1]
			 	end
			end
		end
	else
		for i in 1:length(risk)
			EV_weight-=risk[i][1]
			so=scenario_objs
			if typeof(risk[i]) <: Tuple{Float64,Float64,Dict{Leaf,Float64}}
				for j in 1:length(so)
					so[j]=(so[j][1],so[j][2]-risk[i][3][so[j][3]],so[j][3])
				end
			end
			sort!(so,by=i->i[2],rev=true)
			beta=risk[i][2]
			for scen in so
				if scen[1]>beta
					obj+=scen[2]*risk[i][1]*beta/risk[i][2]
					beta=0
				else
					obj+=scen[2]*risk[i][1]*scen[1]/risk[i][2]
					beta-=scen[1]
				end
			end
		end
	end

	obj+EV_weight*EV
end

function overprint(str)
    print("\e[2K")
    print("\e[1G")
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
