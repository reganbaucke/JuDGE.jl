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
    if !same_shutdowns_at_each_node(sub_problems)
        message *= "Every subproblem must have the same shutdown variables"
        return error(message)
    end
    if !objective_zero(sub_problems)
        message *= "The subproblem objective should be set by @sp_objective(JuMP.Model,AffExpr) not @objective(...)"
        return error(message)
    end
    if !sp_objective_defined(sub_problems)
        message *= "The subproblem objective @sp_objective(JuMP.Model,AffExpr) has not been set to an AffExpr"
        return error(message)
    end

    check_costs(sub_problems)
    check_constraints(sub_problems) #must be at the end of this subroutine!!
end

#### Expansion variables at every node have to be the "same"
function same_expansions_at_each_node(subproblems)
    expansion_structures = collect(map(get_structure_of_expansions, values(subproblems)))
    all(map(x -> x == expansion_structures[1], expansion_structures))
end

function same_shutdowns_at_each_node(subproblems)
    shutdown_structures = collect(map(get_structure_of_shutdowns, values(subproblems)))
    all(map(x -> x == shutdown_structures[1], shutdown_structures))
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
        elseif typeof(sp.ext[:objective_expr])!=AffExpr && typeof(sp.ext[:objective_expr])!=VariableRef
            return false
        end
    end
    true
end

function check_costs(subproblems)
    for (n,sp) in subproblems
        check_sp_costs(sp)
    end
end

function check_sp_costs(model)
    if haskey(model.ext, :expansioncosts) && typeof(model.ext[:expansioncosts])!=AffExpr
        error("@expansioncosts must be provided a linear expression (AffExpr)")
    elseif haskey(model.ext, :maintenancecosts) && typeof(model.ext[:maintenancecosts])!=AffExpr
        error("@maintenancecosts must be provided a linear expression (AffExpr)")
    end
    nothing
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

function get_structure_of_shutdowns(subproblem)
    shutkeys = get_shutdown_keys(subproblem)

    function thinwrapper(var)
        if typeof(var) <: AbstractArray
            collect(eachindex(var))
        else
            ()
        end
    end

    return Dict(i => thinwrapper(subproblem.obj_dict[i]) for i in shutkeys)
end

function get_expansion_keys(model)
    filter(keys(model.obj_dict)) do key
        for (exp,var) in model.ext[:expansions]
            if var === model.obj_dict[key] && !model.ext[:options][exp][1]
                return true
            end
        end
        return false
    end
end

function get_shutdown_keys(model)
    filter(keys(model.obj_dict)) do key
        for (exp,var) in model.ext[:expansions]
            if var === model.obj_dict[key] && model.ext[:options][exp][1]
                return true
            end
        end
        return false
    end
end

function check_constraints(subproblems)
    for (n,sp) in subproblems
        check_sp_constraints(sp)
    end
end

function check_sp_constraints(model)
    loct=list_of_constraint_types(model)
    expansion_keys=get_expansion_keys(model)
    shutdown_keys=get_shutdown_keys(model)

    all_expansion_variables=Array{VariableRef,1}()
    all_shutdown_variables=Array{VariableRef,1}()

    for exp_key in expansion_keys
         var=model[exp_key]
         if typeof(var)==VariableRef
             push!(all_expansion_variables,var)
         elseif typeof(var) <: AbstractArray
             for i in eachindex(var)
                 push!(all_expansion_variables,var[i])
             end
         end
    end

    for shut_key in shutdown_keys
        var=model[shut_key]
        if typeof(var)==VariableRef
            push!(all_shutdown_variables,var)
        elseif typeof(var) <: AbstractArray
            for i in eachindex(var)
                push!(all_shutdown_variables,var[i])
            end
        end
    end

    for ct in loct
        for con in all_constraints(model,ct[1],ct[2])
            con_obj=JuMP.constraint_object(con)
            if typeof(con_obj.func)==VariableRef
                if typeof(con_obj.set)==MathOptInterface.GreaterThan{Float64}
                    if con_obj.func in all_shutdown_variables
                        error("JuDGE Specification Error: Positive coefficient for shutdown variable on LHS of >= constraint")
                    end
                elseif typeof(con_obj.set)==MathOptInterface.LessThan{Float64}
                    if con_obj.func in all_expansion_variables
                        error("JuDGE Specification Error: Positive coefficient for expansion variable on LHS of <= constraint")
                    end
                elseif typeof(con_obj.set)==MathOptInterface.EqualTo{Float64}
                    if con_obj.func in all_expansion_variables || con_obj.func in all_shutdown_variables
                        error("JuDGE Specification Error: Expansion or shutdown variable in == constraint")
                    end
                end
            elseif typeof(con_obj.func) <: GenericAffExpr
                if typeof(con_obj.set)==MathOptInterface.GreaterThan{Float64}
                    for (v,c) in con_obj.func.terms
                        if c>0.0 && v in all_shutdown_variables
                            error("JuDGE Specification Error: Positive coefficient for shutdown variable on LHS of >= constraint")
                        elseif c<0.0 && v in all_expansion_variables
                            error("JuDGE Specification Error: Negative coefficient for expansion variable on LHS of >= constraint")
                        end
                    end
                elseif typeof(con_obj.set)==MathOptInterface.LessThan{Float64}
                    for (v,c) in con_obj.func.terms
                        if c<0.0 && v in all_shutdown_variables
                            error("JuDGE Specification Error: Negative coefficient for shutdown variable on LHS of <= constraint")
                        elseif c>0.0 && v in all_expansion_variables
                            error("JuDGE Specification Error: Positive coefficient for expansion variable on LHS of <= constraint")
                        end
                    end
                elseif typeof(con_obj.set)==MathOptInterface.EqualTo{Float64}
                    for (v,c) in con_obj.func.terms
                        if v in all_shutdown_variables || v in all_expansion_variables
                            error("JuDGE Specification Error: Expansion or shutdown variable in == constraint")
                        end
                    end
                end
            elseif typeof(con_obj.func) <: GenericQuadExpr
                LHS=QuadExpr(AffExpr(0.0),OrderedDict{UnorderedPair{VariableRef},Float64}())
                for (v,c) in con_obj.func.terms
                    if v.a in all_shutdown_variables || v.b in all_shutdown_variables || v.a in all_expansion_variables || v.b in all_expansion_variables
                        error("JuDGE Specification Error: Expansion or shutdown variable in quadratic term")
                    end
                end

                if typeof(con_obj.set)==MathOptInterface.GreaterThan{Float64}
                    for (v,c) in con_obj.func.aff.terms
                        if c>0.0 && v in all_shutdown_variables
                            error("JuDGE Specification Error: Positive coefficient for shutdown variable on LHS of >= constraint")
                        elseif c<0.0 && v in all_expansion_variables
                            error("JuDGE Specification Error: Negative coefficient for expansion variable on LHS of >= constraint")
                        end
                    end
                elseif typeof(con_obj.set)==MathOptInterface.LessThan{Float64}
                    for (v,c) in con_obj.func.aff.terms
                        if c<0.0 && v in all_shutdown_variables
                            error("JuDGE Specification Error: Negative coefficient for shutdown variable on LHS of <= constraint")
                        elseif c>0.0 && v in all_expansion_variables
                            error("JuDGE Specification Error: Positive coefficient for expansion variable on LHS of <= constraint")
                        end
                    end
                elseif typeof(con_obj.set)==MathOptInterface.EqualTo{Float64}
                    for (v,c) in con_obj.func.aff.terms
                        if v in all_shutdown_variables || v in all_expansion_variables
                            error("JuDGE Specification Error: Expansion or shutdown variable in == constraint")
                        end
                    end
                end
            elseif typeof(con_obj.func) == Array{GenericAffExpr{Float64,VariableRef},1}
                for aff in con_obj.func
                    for (v,c) in aff.terms
                        if v in all_shutdown_variables || v in all_expansion_variables
                            error("JuDGE Specification Error: Expansion or shutdown variable in indicator constraint")
                        end
                    end
                end
            else
                warn("JuDGE Specification Error: Unable to verify constraint of type "*string(typeof(con_obj.func)))
            end
        end
    end
    nothing
end
