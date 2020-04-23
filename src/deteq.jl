# Builds a deterministic equivalent by joining the subproblems into
# a single MIP, then adding the constraints that there can only be
# a single investment along any path in the tree.
struct DetEqModel
   problem::JuMP.Model

   function DetEqModel(tree, probability_function, sub_problem_builder, solver)
      println("Establishing deterministic equivalent model for tree: " * string(tree))
      probabilities = probability_function(tree)
      sub_problems = Dict(i => sub_problem_builder(i) for i in collect(tree))
      JuDGE.scale_objectives(sub_problems,probabilities)
      print("Checking sub-problem format...")
      JuDGE.check_specification_is_legal(sub_problems)
      println("Passed")
      print("Building deterministic equivalent problem...")
      problem = build_deteq(sub_problems, tree, probabilities, solver)
      println("Complete")
      return new(problem)
   end
end


function build_deteq(sub_problems, tree::T where T <: AbstractTree, probabilities, solver)
    model = JuMP.Model(solver)

    @objective(model,Min,0)

    history=JuDGE.history(tree)

    model.ext[:vars]=Dict()
    constant=0.0
    for node in keys(sub_problems)
        sp=sub_problems[node]

        objfn=objective_function(sp)
        constant+=objfn.constant

        model.ext[:vars][node] = Dict()
        for variable in all_variables(sp)
            model.ext[:vars][node][string(variable)] = (variable,JuDGE.copy_variable!(model, variable))
            if variable in keys(objfn.terms)
                set_objective_coefficient(model,model.ext[:vars][node][string(variable)][2],objfn.terms[variable])
            end
        end

        for con in all_constraints(sp)
            LHS=AffExpr(0.0)
            for (vv,v) in model.ext[:vars][node]
                try
                    add_to_expression!(LHS,normalized_coefficient(con,v[1]),v[2])
                catch
                    add_to_expression!(LHS,1,v[2])
                end
            end
            set=JuMP.constraint_object(con).set
            if typeof(set)==MathOptInterface.GreaterThan{Float64}
                @constraint(model,LHS>=set.lower)
            elseif typeof(set)==MathOptInterface.LessThan{Float64}
                @constraint(model,LHS<=set.upper)
            elseif typeof(set)==MathOptInterface.EqualTo{Float64}
                @constraint(model,LHS==set.value)
            else
                error("Unsupported constraint type found: "*string(typeof(set)))
            end
        end
    end

    for node in keys(sub_problems)
        sp=sub_problems[node]
        for variable in all_variables(sp)
            for exps in sp.ext[:expansions]
                if isa(exps,VariableRef)
                    if string(exps)==string(variable)
                        model.ext[:vars][node][string(variable)*"_master"] =JuDGE.copy_variable!(model, variable)
                        set_objective_coefficient(model,model.ext[:vars][node][string(variable)*"_master"],probabilities(node)*JuDGE.coef(sp.ext[:expansioncosts],variable))
                    end
                elseif typeof(exps) <: AbstractArray
                    for expvar in exps
                        if string(expvar)==string(variable)
                            model.ext[:vars][node][string(variable)*"_master"] =JuDGE.copy_variable!(model, variable)
                            set_objective_coefficient(model,model.ext[:vars][node][string(variable)*"_master"],probabilities(node)*JuDGE.coef(sp.ext[:expansioncosts],variable))
                        end
                    end
                end
            end
        end
    end

    for node in keys(sub_problems)
        sp=sub_problems[node]
        for variable in all_variables(sp)
            for exps in sp.ext[:expansions]
                if isa(exps,VariableRef)
                    if string(exps)==string(variable)
                        @constraint(model,model.ext[:vars][node][string(variable)][2]<=sum(model.ext[:vars][n][string(variable)*"_master"] for n in history(node)))
                        if typeof(node)==Leaf
                            @constraint(model,sum(model.ext[:vars][n][string(variable)*"_master"] for n in history(node))<=1)
                        end
                    end
                elseif typeof(exps) <: AbstractArray
                    for expvar in exps
                        if string(expvar)==string(variable)
                            @constraint(model,model.ext[:vars][node][string(variable)][2]<=sum(model.ext[:vars][n][string(variable)*"_master"] for n in history(node)))
                            if typeof(node)==Leaf
                                @constraint(model,sum(model.ext[:vars][n][string(variable)*"_master"] for n in history(node))<=1)
                            end
                        end
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
