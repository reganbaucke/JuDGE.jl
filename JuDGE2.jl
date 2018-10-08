module JuDGE2

push!(LOAD_PATH, ".")

using JuDGETree
using JuMP
# using Gurobi


mutable struct JuDGEModel
    tree::Tree
    master::JuMP.Model
    subprob::Dict{Node,JuMP.Model}
    duals::Dict{Node,Any}
    mastervar::Dict{Node,Dict{Symbol,Any}}
    mastercon::Dict{Node,Dict{Symbol,Any}}
    buildexpansionvariables
    buildsubs
    expansioncosts
    function JuDGEModel(tree::Tree)
        this = new()
        this.tree = tree
        return this
    end
end

function JuDGEexpansions!(f!,jmodel::JuDGEModel)
    jmodel.buildexpansionvariables = f!
end

function JuDGEsubproblem!(f!,jmodel::JuDGEModel)
    jmodel.buildsubs = f!
end

function JuDGEexpansioncosts!(f!,jmodel::JuDGEModel)
    jmodel.expansioncosts = f!
end
 
function buildsubproblems(jmodel::JuDGEModel)
    jmodel.subprob=Dict{Node,JuMP.Model}()
    for n in jmodel.tree.nodes
        sp = JuMP.Model()
        jmodel.buildexpansionvariables(sp)
        jmodel.buildsubs(sp,n)
        jmodel.subprob[n] = sp
    end
end

function buildmaster(jmodel::JuDGEModel)
    # create the jump model
    jmodel.master = JuMP.Model()

    # initialize the dicts
    mastervar = Dict{Node,Dict{Symbol,Any}}()
    mastercon = Dict{Node,Dict{Symbol,Any}}()
    for n in jmodel.tree.nodes
        mastervar[n] = Dict{Symbol,Any}()
        mastercon[n] = Dict{Symbol,Any}()
    end

    sp = jmodel.subprob[jmodel.tree.root]
    for (key,value) in filter((key,value) -> value == :expansion,sp.ext)
        println(key)
        if isa(sp.objDict[key], JuMP.Variable)
            println("single")
            for n in jmodel.tree.nodes
                mastervar[n][key] = @variable(jmodel.master)
            end
        elseif isa(sp.objDict[key], AbstractArray)
            println("normal array")
            for n in jmodel.tree.nodes
                mastervar[n][key] = @variable(jmodel.master)
            end
        elseif isa(sp.objDict[key], JuMP.JuMPArray)
            println("jump array")
            for n in jmodel.tree.nodes
                mastervar[n][key] = @variable(jmodel.master,[])
            end
        end
    end
end

function JuDGEbuild!(jmodel::JuDGEModel)
    buildsubproblems(jmodel)
    buildmaster(jmodel)
end

function solve(jmodel::JuDGEModel,iter::Int64)
    # perform an iteration
    for i = 1:iter
        iteration(jmodel)
        println("$i/$iter")
    end
end

function iteration(jmodel::JuDGEModel)
    solved = JuMP.solve(jmodel.master)
    println("Master problem objective function")
    println(jmodel.master.objVal)
    for n in jmodel.tree.nodes
        if solved == :Optimal
            jmodel.updateduals(n)
        end
        jmodel.updateobjective(n)
        JuMP.solve(jmodel.subprob[n])
        addcolumn(jmodel, jmodel.buildcolumn(n))
    end
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

# macro dual(pi, indices...)
#     tmp = String(pi) * " = makeDual( "
#     for set in indices[1:end-1]
#         if typeof(set) == Symbol
#             tmp = tmp * String(set) * ", "
#         elseif typeof(set) == Expr
#             tmp = tmp * repr(set)[3:end-1] * ", "
#         end
#     end
#     if typeof(indices[end]) == Symbol
#         tmp = tmp * String(indices[end]) *")"
#     elseif typeof(indices[end]) == Expr
#         tmp = tmp * repr(indices[end])[3:end-1] *")"
#     end
#     tmp2 = String(pi) * " = (:" * String(pi) * " , " * String(pi) * ")"
#     final = quote
#         $(esc(parse(tmp)))
#         $(esc(parse(tmp2)))
#     end
#     return final
# end

# macro dual(pi)
#     tmp = String(pi) * " = makeDual() "
#     tmp2 = String(pi) * " = (:" * String(pi) * " , " * String(pi) * ")"
#     final = quote
#         $(esc(parse(tmp)))
#         $(esc(parse(tmp2)))
#     end
#     return final
# end

macro expansion(sp,x::Expr)
    tmp = "@variable(" * String(sp) * "," * String(repr(x))[3:end-1] * ", category=:Bin)"
    tmp2 = "sp.ext[:" * String(x.args[1]) * "] = :expansion"
    final = quote
        $(esc(parse(tmp)))
        $(esc(parse(tmp2)))
    end
    return final
end

macro expansion(sp,x::Symbol)
    tmp = "@variable(" * String(sp) * "," * String(x) * ", category=:Bin)"
    tmp2 = "sp.ext[:" * String(x) * "] = :expansion"
    final = quote
        $(esc(parse(tmp2)))
        $(esc(parse(tmp)))
    end
    return final
end

macro bullshit(sp,x::Expr)
    # x =:($x)
    println(x.head)
    println(length(x.args))
    println(x.args)
    tmp = "@variable(" * String(sp) * "," * String(repr(x))[3:end-1][1] * ", category=:Bin)"
    final = quote
        # $(esc(parse(tmp2)))
        $(esc(parse(tmp)))
    end
    return final
end

# end module
export 
JuDGEsubproblem!, JuDGEModel, solve, JuDGEexpansions!, @expansion, JuDGEbuild!

end