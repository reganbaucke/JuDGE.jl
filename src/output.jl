function JuDGE_value(jmodel::JuDGEModel, node::AbstractTree, var::Symbol, index::Int64)
   if termination_status(jmodel.master_problem) != MathOptInterface.OPTIMAL
      error("JuDGE model not solved")
   end
   return JuMP.value(jmodel.sub_problems[node][var][index])
end

function JuDGE_value(jmodel::JuDGEModel, node::AbstractTree, var::Symbol)
   if termination_status(jmodel.master_problem) != MathOptInterface.OPTIMAL
      error("JuDGE model not solved")
   end
   return JuMP.value(jmodel.sub_problems[node][var])
end

function print_expansions(jmodel::JuDGEModel;node=jmodel.tree::AbstractTree,onlynonzero::Bool=true)
    if termination_status(jmodel.master_problem) != MathOptInterface.OPTIMAL
        error("You need to first solve the decomposed model.")
    end
    # this is how you access the value of the binary expansions in the master
    for x in keys(jmodel.master_problem.ext[:expansions][node])
        var = jmodel.master_problem.ext[:expansions][node][x]
        # print value of variable (if non-zero)
        if isa(var,Array)
            for key in keys(var)
                if !onlynonzero || value(var[key...])>0
                    println(node.name* "_" * string(x) * "[" * string(key)* "]" * ": " * string(value(var[key...])))
                end
            end
        else
            println(node.name * "_" * string(x) *": " * string(value(var)))
        end
    end

    if  typeof(node)==Tree && !(length(node.children) == 0)
        for i in 1:length(node.children)
            print_expansions(jmodel,node=node.children[i],onlynonzero=onlynonzero)
        end
    end
end
