###
# asserts for constructor of JuDGEModel
###

# Definition of an important set of constraints
function NormalConstraints()
    [
        (GenericAffExpr{Float64,VariableRef}, MathOptInterface.LessThan{Float64}),
        (GenericAffExpr{Float64,VariableRef}, MathOptInterface.GreaterThan{Float64}),
        (GenericAffExpr{Float64,VariableRef}, MathOptInterface.EqualTo{Float64}),
        (VariableRef, MathOptInterface.LessThan{Float64}),
        (VariableRef, MathOptInterface.GreaterThan{Float64}),
        (VariableRef, MathOptInterface.EqualTo{Float64}),
    ]
end

function check_specification_is_legal(sub_problems)
    message = "JuDGE Specification Error: "
    if !same_expansions_at_each_node(sub_problems)
        message *= "Every subproblem must have the same expansion variables"
        return error(message)
    end
    if !all(map(expansions_only_in_expansion_constraints, values(sub_problems)))
        message *= "Expansion variables can only appear in expansion constraints"
        return error(message)
    end
    # if !same_expansion_constraints_at_each_node(sub_problems)
    #     message *= "Every subproblem must have the same expansion constraints (but not neccesarily the same coefficients)"
    #     return error(message)
    # end
    if !expansion_constraints_are_less_than(sub_problems)
        message *= "Expansion constraints must be \'less than or equal to\' constraints"
        return error(message)
    end
    #if !expansion_variables_same_coeffs_in_expansion_constraints(sub_problems)
    #    message *= "Expansion variables in expansion constraints must be have the same coefficients between all subproblems"
    #    return error(message)
    #end
    if !coeffs_are_positive(sub_problems)
        message *= "Coefficients of expansion variables in the expansion constraints must be non-negative and on the RHS (alternatively non-positive and on the LHS)"
        return error(message)
    end
    # if !same_constant_terms(sub_problems)
    #     message *= "The constant term in each expansion constraint must be the same value between all subproblems"
    #     return error(message)
    # end
    if !objective_zero(sub_problems)
        message *= "The subproblem objective should be set by @sp_objective(JuMP.Model,AffExpr) not @objective(...)"
        return error(message)
    end
    if !sp_objective_defined(sub_problems)
        message *= "The subproblem objective @sp_objective(JuMP.Model,AffExpr) has not been set"
        return error(message)
    end
    return nothing
end

#### Subproblems all must have the same sense
function same_sense_at_each_node(subproblems)
    sense_of_subproblems = collect(map(objective_sense, values(subproblems)))
    all(map(x -> x == MathOptInterface.OptimizationSense(0), sense_of_subproblems))
end

#### Expansion variables at every node have to be the "same"
function same_expansions_at_each_node(subproblems)
    expansion_structures = collect(map(get_structure_of_expansions, values(subproblems)))
    all(map(x -> x == expansion_structures[1], expansion_structures))
end

#### Expansion constraints has to have have to be the "same"
function same_expansion_constraints_at_each_node(subproblems)
    ec_of_subproblems =
        collect(map(get_structure_of_expansion_constraints, values(subproblems)))
    all(map(x -> x == ec_of_subproblems[1], ec_of_subproblems))
end

#### Expansion constraints have to be less than or equal to
function expansion_constraints_are_less_than(subproblems)
    less_than = collect(map(are_less_than, values(subproblems)))
    f(x) = all(map(all, values(map(all, x))))
    all(map(f, less_than))
end

#### Expansion variables have to to have same coeffs in the expansion constraints
function expansion_variables_same_coeffs_in_expansion_constraints(subproblems)
    coeffs = collect(map(get_constraints_coeffs_for_expansions, values(subproblems)))
    all(map(x -> x == coeffs[1], coeffs))
end

#### initial capacity is the same at every node
function same_constant_terms(subproblems)
    constant_terms = collect(map(get_constant_terms, values(subproblems)))
    all(map(x -> x == constant_terms[1], constant_terms))
end

function objective_zero(subproblems)
    objective_functions = collect(map(objective_function, values(subproblems)))
    for objfunc in objective_functions
        if length(objfunc.terms)!=0
            return false
        end
    end
    true
end

function sp_objective_defined(subproblems)
    for (n,sp) in subproblems
        if !haskey(sp.ext, :objective)
            return false
        end
    end
    true
end

#### coeffs of expansion variables are non-negative in the expansion constraints
function coeffs_are_positive(subproblems)
    nodes = collect(map(coeffs_of_expansion_variables_are_positive, values(subproblems)))

    # complicated way of checking that every coeff is positive
    # for a single subproblem, each coeff is positive
    f(x) = all(map(all, values(map(all, x))))

    # check that this is true for all subs
    all(map(f, nodes))
end

function expansion_variables_in_expansion_constraints_only(subproblems)
    nothing
end

function is_less_than(
    ::ConstraintRef{
        Model,
        MathOptInterface.ConstraintIndex{
            MathOptInterface.ScalarAffineFunction{Float64},
            MathOptInterface.LessThan{Float64},
        },
        ScalarShape,
    },
)
    true
end
function is_less_than(other)
    false
end

function get_structure_of_expansions(subproblem)
    expkeys = get_expansion_keys(subproblem)

    function thinwrapper(var)
        if typeof(var) <: AbstractArray
            collect(eachindex(var))
        else
            ()
        end
    end

    return Dict(i => thinwrapper(subproblem.obj_dict[i]) for i in expkeys)
end

function get_constraints_coeffs_for_expansions(subproblem)
    expkeys = get_expansion_keys(subproblem)
    expconkeys = get_expansion_constraint_keys(subproblem)

    function thinwrapper(var, con)
        if typeof(var) <: AbstractArray
            if typeof(con) <: AbstractArray
                holder = [(i, j) for i in eachindex(var), j in eachindex(con)]
                map(holder) do (i, j)
                    normalized_coefficient(con[j], var[i])
                end
            else
                holder = [i for i in eachindex(var)]
                map(holder) do i
                    normalized_coefficient(con, var[i])
                end
            end
        else
            if typeof(con) <: AbstractArray
                holder = [j for jj in eachindex(con)]
                map(holder) do j
                    normalized_coefficient(con, var[j])
                end
            else
                normalized_coefficient(con, var)
            end
        end
    end

    return Dict(
        (i, j) => thinwrapper(subproblem.obj_dict[i], subproblem.obj_dict[j])
        for i in expkeys, j in expconkeys
    )
end

function coeffs_of_expansion_variables_are_positive(subproblem)
    expkeys = get_expansion_keys(subproblem)
    expconkeys = get_expansion_constraint_keys(subproblem)

    function thinwrapper(var, con)
        if typeof(var) <: AbstractArray
            if typeof(con) <: AbstractArray
                holder = [(i, j) for i in eachindex(var), j in eachindex(con)]
                map(holder) do (i, j)
                    normalized_coefficient(con[j], var[i]) <= 0.0
                end
            else
                holder = [i for i in eachindex(var)]
                map(holder) do i
                    normalized_coefficient(con, var[i]) <= 0.0
                end
            end
        else
            if typeof(con) <: AbstractArray
                holder = [j for j in eachindex(con)]
                map(holder) do j
                    normalized_coefficient(con, var[j]) <= 0.0
                end
            else
                normalized_coefficient(con, var) <= 0.0
            end
        end
    end

    return Dict(
        (i, j) => thinwrapper(subproblem.obj_dict[i], subproblem.obj_dict[j])
        for i in expkeys, j in expconkeys
    )
end

function get_structure_of_expansion_constraints(subproblem)
    expconkeys = get_expansion_constraint_keys(subproblem)

    function thinwrapper(var)
        if typeof(var) <: AbstractArray
            collect(eachindex(var))
        else
            ()
        end
    end

    return Dict(i => thinwrapper(subproblem.obj_dict[i]) for i in expconkeys)
end

function are_less_than(subproblem)
    expconkeys = get_expansion_constraint_keys(subproblem)

    function thinwrapper(con)
        if typeof(con) <: AbstractArray
            map(is_less_than, con)
        else
            is_less_than(con)
        end
    end

    return Dict(i => thinwrapper(subproblem.obj_dict[i]) for i in expconkeys)
end

function get_constant_terms(subproblem)
    expconkeys = get_expansion_constraint_keys(subproblem)

    function thinwrapper(con)
        if typeof(con) <: AbstractArray
            map(normalized_rhs, con)
        else
            normalized_rhs(con)
        end
    end
    return Dict(i => thinwrapper(subproblem.obj_dict[i]) for i in expconkeys)
end

function expansions_only_in_expansion_constraints(subproblem)
    expkeys = get_expansion_keys(subproblem)
    expconkeys = get_expansion_constraint_keys(subproblem)

    for con in setdiff(
        Set(all_constraints(subproblem)),
        Set(semicollect(Dict(i => subproblem.obj_dict[i] for i in expconkeys))),
    )
        expr = constraint_object(con).func
        if typeof(expr) == GenericAffExpr{Float64,VariableRef}
            variables_in_con = linear_terms(constraint_object(con).func)
            variables_in_con = map(x -> x[2], variables_in_con)
        elseif typeof(expr) == VariableRef
            variables_in_con = [expr]
        else
            error("Cannot check this constraint type")
        end

        for var in semicollect(Dict(i => subproblem.obj_dict[i] for i in expkeys))
            if var in variables_in_con
                return false
            end
        end
    end
    return true
end

function wrapper(input::T where {T<:AbstractArray})
    input
end

function wrapper(input)
    [input]
end

function semicollect(dict::Dict)
    out = []
    for key in keys(dict)
        for el in wrapper(dict[key])
            push!(out, el)
        end
    end
    out
end

function JuMP.all_constraints(model::JuMP.Model)
    types = list_of_constraint_types(model)
    out = []
    for type in NormalConstraints()
        for con in all_constraints(model, type...)
            push!(out, con)
        end
    end
    out
end

function get_expansion_keys(model)
    filter(keys(model.obj_dict)) do key
        for (exp,var) in model.ext[:expansions]
            if var === model.obj_dict[key] && !model.ext[:forced][exp]
                return true
            end
        end
        return false
    end
end

function get_expansion_constraint_keys(model)
    filter(keys(model.obj_dict)) do key
        for exp in model.ext[:expansionconstraints]
            if exp === model.obj_dict[key]
                return true
            end
        end
        return false
    end
end
