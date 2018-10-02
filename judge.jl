module JuDGE

push!(LOAD_PATH, ".")

using JuDGETree
using JuMP
using Gurobi


mutable struct JuDGEModel
    tree::Tree
    master::JuMP.Model
    subprob::Dict{Node,JuMP.Model}
    duals::Dict{Node,Any}
    updateduals
    buildcolumn
    updateobjective
    function JuDGEModel(tree::Tree)
        this = new()
        this.tree = tree
        return this
    end
end

mutable struct Dual
    indexset::Array{Any,1}
    value::Array{Float64,1}
end


function JuDGEduals!(f,jmodel::JuDGEModel)
    jmodel.duals = Dict{Node,Any}()
    for n in jmodel.tree.nodes
        jmodel.duals[n] = Dict{Symbol,Dual}()
        out = f(n)
        for dual in out
            jmodel.duals[n][dual[1]] = dual[2]
        end
    end
end

function JuDGEsubproblems!(f,jmodel::JuDGEModel)
    jmodel.subprob = Dict{Node,JuMP.Model}()
    for n in jmodel.tree.nodes
        jmodel.subprob[n] = f(n)
    end

end

function JuDGEobjective!(f,jmodel::JuDGEModel)
    jmodel.updateobjective = f
end

function JuDGEupdateduals!(f,jmodel::JuDGEModel)
    jmodel.updateduals = f
end

function JuDGEbuildcolumn!(f,jmodel::JuDGEModel)
    jmodel.buildcolumn = f
end

function JuDGEmaster!(f,jmodel::JuDGEModel)
    jmodel.master = f()
end

function addcolumn(jmodel,column)
    Variable(jmodel.master,column...)
end

function solve(jmodel::JuDGEModel)
    # perform an iteration
    for i = 1:5
        iteration(jmodel)
    end
end

function iteration(jmodel::JuDGEModel)
    solved = JuMP.solve(jmodel.master)
    for n in jmodel.tree.nodes
        if solved == :Optimal
            jmodel.updateduals(n)
        end
        jmodel.updateobjective(n)
        JuMP.solve(jmodel.subprob[n])
        addcolumn(jmodel, jmodel.buildcolumn(n))
    end
end

function makeDual(iterables...)
    # because the dictionary breaks, here we have to see if iterable is unique and use an array
    # indexset = Array{Any,1}()
    indexset = [];
    value = Array{Float64,1}()
    for i in Iterators.product(iterables...)
        push!(indexset,i)
        push!(value,0.0)
    end
    Dual(indexset,value)
end

function makeDual()
    # because the dictionary breaks, here we have to see if iterable is unique and use an array
    # indexset = Array{Any,1}()
    indexset = [()];
    value = [0.0]
    return Dual(indexset,value)
end

function Base.getindex(pi::Dual, i...)
    pi.value[findfirst(pi.indexset,i)]
end

function Base.setindex!(pi::Dual, value ,i...)
    pi.value[findfirst(pi.indexset,i)] = value
end

function JuDGEModel(f,tree::Tree)
    master = JuMP.Model(solver=GurobiSolver())
    subprob = Dict{Node{Tree},JuMP.Model}()
    duals = Dict{Node{Tree},Any}()
    for n in tree.nodes
        # sp = JuMP.Model(solver=GurobiSolver())
        # subprob[n] = f(n)
        subprob[n] = f(n)
    end
    JuDGEModel(tree,master,subprob,duals,nothing,nothing)
end

macro dual(pi, indices...)
    tmp = String(pi) * " = makeDual( "
    for set in indices[1:end-1]
        if typeof(set) == Symbol
            tmp = tmp * String(set) * ", "
        elseif typeof(set) == Expr
            tmp = tmp * repr(set)[3:end-1] * ", "
        end
    end
    if typeof(indices[end]) == Symbol
        tmp = tmp * String(indices[end]) *")"
    elseif typeof(indices[end]) == Expr
        tmp = tmp * repr(indices[end])[3:end-1] *")"
    end
    tmp2 = String(pi) * " = (:" * String(pi) * " , " * String(pi) * ")"
    final = quote
        $(esc(parse(tmp)))
        $(esc(parse(tmp2)))
    end
    return final
end

macro dual(pi)
    tmp = String(pi) * " = makeDual() "
    tmp2 = String(pi) * " = (:" * String(pi) * " , " * String(pi) * ")"
    final = quote
        $(esc(parse(tmp)))
        $(esc(parse(tmp2)))
    end
    return final
end

# end module
export 
JuDGEbuildcolumn!, JuDGEduals!, JuDGEmaster!, JuDGEobjective!, 
JuDGEsubproblems!, JuDGEModel, JuDGEupdateduals!, @dual, makeDual, solve

end