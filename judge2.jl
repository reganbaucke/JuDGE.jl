using JuMP

mutable struct JuDGEModel
    tree::Tree
    master::JuMP.Model
    subprob::Dict{Node,JuMP.Model}
end

function subproblem(f::Function, model::JuDGEModel)
    for n in model.tree.nodes
        model.subprob[n] = f(n)
    end
end
