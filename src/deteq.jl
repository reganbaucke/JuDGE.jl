# Builds a deterministic equivalent by joining the subproblems into
# a single MIP, then adding the constraints that there can only be
# a single investment along any path in the tree.
using DataStructures

struct DetEqModel
   problem::JuMP.Model
   function DetEqModel(tree, probability_function, sub_problem_builder, solver; discount_factor=1.0)
      println("Establishing deterministic equivalent model for tree: " * string(tree))
      probabilities = probability_function(tree)
      sub_problems = Dict(i => sub_problem_builder(i) for i in collect(tree))
      JuDGE.scale_objectives(tree,sub_problems,probabilities,discount_factor)
      print("Checking sub-problem format...")
      JuDGE.check_specification_is_legal(sub_problems)
      println("Passed")
      print("Building deterministic equivalent problem...")
      problem = build_deteq(sub_problems, tree, probabilities, solver, discount_factor)
      println("Complete")
      return new(problem)
   end
end

function build_deteq(sub_problems, tree::T where T <: AbstractTree, probabilities, solver, discount_factor::Float64)
    model = JuMP.Model(solver)

    @objective(model,Min,0)

    history=JuDGE.history(tree)
    depth=JuDGE.depth(tree)

    model.ext[:vars]=Dict()
    constant=0.0
    for (node,sp) in sub_problems
        objfn=objective_function(sp)
        constant+=objfn.constant

        model.ext[:vars][node] = Dict()
        for variable in all_variables(sp)
            model.ext[:vars][node][variable] = JuDGE.copy_variable!(model, variable)
            if variable in keys(objfn.terms)
                set_objective_coefficient(model,model.ext[:vars][node][variable],objfn.terms[variable])
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
                elseif typeof(set)!=MathOptInterface.ZeroOne
                    error("Unsupported constraint type found: "*string(typeof(set)))
                else
                    continue
                end
            end
        end
    end

    for (node,sp) in sub_problems
        df=discount_factor^depth(node)
        for (name,exps) in sp.ext[:expansions]
            if isa(exps,VariableRef)
                variable=sp[name]
                model.ext[:vars][node][string(variable)*"_master"] =JuDGE.copy_variable!(model, variable)
                set_objective_coefficient(model,model.ext[:vars][node][string(variable)*"_master"],df*probabilities(node)*JuDGE.coef(sp.ext[:expansioncosts],variable))
            elseif typeof(exps) <: AbstractArray
                variables=sp[name]
                for index in eachindex(exps)
                    model.ext[:vars][node][string(variables[index])*"_master"] =JuDGE.copy_variable!(model, variables[index])
                    set_objective_coefficient(model,model.ext[:vars][node][string(variables[index])*"_master"],df*probabilities(node)*JuDGE.coef(sp.ext[:expansioncosts],variables[index]))
                end
            end
        end
    end

    for (node,sp) in sub_problems
        for (name,exps) in sp.ext[:expansions]
            if isa(exps,VariableRef)
                @constraint(model,model.ext[:vars][node][exps]<=sum(model.ext[:vars][n][string(exps)*"_master"] for n in history(node)))
                if typeof(node)==Leaf
                    @constraint(model,sum(model.ext[:vars][n][string(exps)*"_master"] for n in history(node))<=1)
                end
            elseif typeof(exps) <: AbstractArray
                for expvar in exps
                    @constraint(model,model.ext[:vars][node][expvar]<=sum(model.ext[:vars][n][string(expvar)*"_master"] for n in history(node)))
                    if typeof(node)==Leaf
                        @constraint(model,sum(model.ext[:vars][n][string(expvar)*"_master"] for n in history(node))<=1)
                    end
                end
            end
        end
    end
    @objective(model,Min,objective_function(model)+constant)
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
