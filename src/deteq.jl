# Builds a deterministic equivalent by joining the subproblems into
# a single MIP, then adding the constraints that there can only be
# a single investment along any path in the tree.
using DataStructures

struct DetEqModel
   problem::JuMP.Model
   function DetEqModel(tree, probability_function, sub_problem_builder, solver; discount_factor=1.0, CVaR=(0.0,1.0), intertemporal=nothing)
      println("")
      println("Establishing deterministic equivalent model for tree: " * string(tree))
      probabilities = probability_function(tree)
      sub_problems = Dict(i => sub_problem_builder(i) for i in collect(tree))
      print("Checking sub-problem format...")
      JuDGE.check_specification_is_legal(sub_problems)
      println("Passed")
      JuDGE.scale_objectives(tree,sub_problems,discount_factor)
      print("Building deterministic equivalent problem...")
      problem = build_deteq(sub_problems, tree, probabilities, solver, discount_factor, CVaR, intertemporal)
      println("Complete")
      return new(problem)
   end
end

function build_deteq(sub_problems, tree::T where T <: AbstractTree, probabilities, solver, discount_factor::Float64, CVaR::Tuple{Float64,Float64}, intertemporal)
    model = JuMP.Model(solver)

    @objective(model,Min,0)

    history_function=JuDGE.history(tree)
    depth_function=JuDGE.depth(tree)

    model.ext[:vars]=Dict()
    model.ext[:master_vars]=Dict{AbstractTree,Dict{Symbol,Any}}()
    model.ext[:master_names]=Dict{AbstractTree,Dict{Symbol,Any}}()

    scen_con=Dict{Leaf,ConstraintRef}()
    scen_var=Dict{Leaf,VariableRef}()

    for leaf in get_leafnodes(tree)
        scen_var[leaf]=@variable(model)
        scen_con[leaf]=@constraint(model,0==scen_var[leaf])
        set_objective_coefficient(model,scen_var[leaf],probabilities(leaf))
    end

    if CVaR[1]>0.0 && CVaR[2]<1.0
        eta=@variable(model)
        for leaf in get_leafnodes(tree)
            v = @variable(model)
            w = @variable(model)
            set_lower_bound(v,0.0)
            set_lower_bound(w,0.0)
            @constraint(model,v>=eta-scen_var[leaf])
            @constraint(model,w>=scen_var[leaf]-eta)
            set_objective_coefficient(model, v, probabilities(leaf)*CVaR[1])
            set_objective_coefficient(model, w, probabilities(leaf)*CVaR[1]/CVaR[2]*(1-CVaR[2]))
        end
    end

    for (node,sp) in sub_problems
        #objfn=objective_function(sp)
        #constant+=objfn.constant
        leafnodes=get_leafnodes(node)
        model.ext[:vars][node] = Dict()
        for variable in all_variables(sp)
            model.ext[:vars][node][variable] = JuDGE.copy_variable!(model, variable)
            if variable==sp.ext[:objective]
                for leaf in leafnodes
                    set_normalized_coefficient(scen_con[leaf],model.ext[:vars][node][variable],1.0)
                end
            end
        end

        loct=list_of_constraint_types(sp)

        for ct in loct
            for con in all_constraints(sp,ct[1],ct[2])
                con_obj=JuMP.constraint_object(con)
                if typeof(con_obj.func)==VariableRef
                    LHS=AffExpr(0.0)
                    add_to_expression!(LHS,1,model.ext[:vars][node][JuMP.constraint_object(con).func])
                elseif typeof(con_obj.func) <: GenericAffExpr
                    LHS=AffExpr(0.0)
                    for (v,c) in con_obj.func.terms
                        add_to_expression!(LHS,c,model.ext[:vars][node][v])
                    end
                elseif typeof(con_obj.func) <: GenericQuadExpr
                    LHS=QuadExpr(AffExpr(0.0),OrderedDict{UnorderedPair{VariableRef},Float64}())
                    for (v,c) in con_obj.func.terms
                        LHS.terms[UnorderedPair{VariableRef}(model.ext[:vars][node][v.a],model.ext[:vars][node][v.b])]=c
                    end
                    for (v,c) in con_obj.func.aff.terms
                        add_to_expression!(LHS.aff,c,model.ext[:vars][node][v])
                    end
                elseif typeof(con_obj.func) == Array{GenericAffExpr{Float64,VariableRef},1}
                    group=Array{AffExpr,1}()
                    for aff in con_obj.func
                        LHS=AffExpr(0.0)
                        for (v,c) in aff.terms
                            add_to_expression!(LHS,c,model.ext[:vars][node][v])
                        end
                        push!(group,LHS)
                    end
                else
                    error("Unsupported constraint type found: "*string(typeof(con_obj.func)))
                end
                set=con_obj.set
                if typeof(set)==MathOptInterface.GreaterThan{Float64}
                    @constraint(model,LHS>=set.lower)
                elseif typeof(set)==MathOptInterface.LessThan{Float64}
                    @constraint(model,LHS<=set.upper)
                elseif typeof(set)==MathOptInterface.EqualTo{Float64}
                    @constraint(model,LHS==set.value)
                elseif typeof(set)==MathOptInterface.SecondOrderCone
                    @constraint(model,group in SecondOrderCone())
                elseif typeof(set)==MathOptInterface.IndicatorSet{MathOptInterface.ACTIVATE_ON_ZERO,MathOptInterface.EqualTo{Float64}}
                    @constraint(model,!collect(keys(group[1].terms))[1] => {group[2]==set.set.value})
                elseif typeof(set)==MathOptInterface.IndicatorSet{MathOptInterface.ACTIVATE_ON_ZERO,MathOptInterface.LessThan{Float64}}
                    @constraint(model,!collect(keys(group[1].terms))[1] => {group[2]<=set.set.value})
                elseif typeof(set)==MathOptInterface.IndicatorSet{MathOptInterface.ACTIVATE_ON_ZERO,MathOptInterface.GreaterThan{Float64}}
                    @constraint(model,!collect(keys(group[1].terms))[1] => {group[2]>=set.set.value})
                elseif typeof(set)==MathOptInterface.IndicatorSet{MathOptInterface.ACTIVATE_ON_ONE,MathOptInterface.EqualTo{Float64}}
                    @constraint(model,collect(keys(group[1].terms))[1] => {group[2]==set.set.value})
                elseif typeof(set)==MathOptInterface.IndicatorSet{MathOptInterface.ACTIVATE_ON_ONE,MathOptInterface.LessThan{Float64}}
                    @constraint(model,collect(keys(group[1].terms))[1] => {group[2]<=set.set.value})
                elseif typeof(set)==MathOptInterface.IndicatorSet{MathOptInterface.ACTIVATE_ON_ONE,MathOptInterface.GreaterThan{Float64}}
                    @constraint(model,collect(keys(group[1].terms))[1] => {group[2]>=set.set.value})
                elseif typeof(set)!=MathOptInterface.ZeroOne
                    error("Unsupported constraint type found: "*string(typeof(set)))
                else
                    continue
                end
            end
        end
    end

    for (node,sp) in sub_problems
        model.ext[:master_vars][node] = Dict{Symbol,Any}()
        model.ext[:master_names][node] = Dict{Symbol,Any}()
        leafnodes=get_leafnodes(node)
        df=discount_factor^depth_function(node)
        for (name,exps) in sp.ext[:expansions]
            if isa(exps,VariableRef)
                variable=sp[name]
                model.ext[:master_vars][node][name] =JuDGE.copy_variable!(model, variable)
                model.ext[:master_names][node][name] = string(variable)
                for leaf in leafnodes
                    set_normalized_coefficient(scen_con[leaf],model.ext[:master_vars][node][name],df*JuDGE.coef(sp.ext[:expansioncosts],variable))
                end
                #set_objective_coefficient(model,model.ext[:vars][node][string(variable)*"_master"],df*probabilities(node)*JuDGE.coef(sp.ext[:expansioncosts],variable))
            elseif typeof(exps) <: AbstractArray
                variables=sp[name]
                model.ext[:master_vars][node][name]=Dict()
                model.ext[:master_names][node][name]=Dict()
                for index in eachindex(exps)
                    model.ext[:master_vars][node][name][index] =JuDGE.copy_variable!(model, variables[index])
                    model.ext[:master_names][node][name][index]=string(variables[index])
                    for leaf in leafnodes
                        set_normalized_coefficient(scen_con[leaf],model.ext[:master_vars][node][name][index],df*JuDGE.coef(sp.ext[:expansioncosts],variables[index]))
                    end
                    #set_objective_coefficient(model,model.ext[:vars][node][string(variables[index])*"_master"],df*probabilities(node)*JuDGE.coef(sp.ext[:expansioncosts],variables[index]))
                end
            end
        end
    end

    for (node,sp) in sub_problems
        for (name,exps) in sp.ext[:expansions]
            if isa(exps,VariableRef)
                if sp.ext[:forced][name]
                    @constraint(model,model.ext[:vars][node][exps]==sum(model.ext[:master_vars][n][name] for n in history_function(node)))
                else
                    @constraint(model,model.ext[:vars][node][exps]<=sum(model.ext[:master_vars][n][name] for n in history_function(node)))
                end
                if typeof(node)==Leaf
                    @constraint(model,sum(model.ext[:master_vars][n][name] for n in history_function(node))<=1)
                end
            elseif typeof(exps) <: AbstractArray
                for index in eachindex(exps)
                    if sp.ext[:forced][name]
                        @constraint(model,model.ext[:vars][node][exps[index]]==sum(model.ext[:master_vars][n][name][index] for n in history_function(node)))
                    else
                        @constraint(model,model.ext[:vars][node][exps[index]]<=sum(model.ext[:master_vars][n][name][index] for n in history_function(node)))
                    end
                    if typeof(node)==Leaf
                        @constraint(model,sum(model.ext[:master_vars][n][name][index] for n in history_function(node))<=1)
                    end
                end
            end
        end
    end

    if typeof(intertemporal) <: Function
        map(Main.eval,unpack_expansions(model.ext[:master_vars])) #bring expansion variables into global scope
        intertemporal(model,tree)
        map(Main.eval,clear_expansions(model.ext[:master_vars]))
    end
    return model
end

# Function called by the user to solve the deterministic equivalent
function solve(deteq::DetEqModel)
    print("Solving deterministic equivalent formulation...")
    optimize!(deteq.problem)
    if termination_status(deteq.problem) == MathOptInterface.OPTIMAL
        println("Solved.")
    else
        println("Not solved: " * string(termination_status(deteq.problem)))
    end
end
