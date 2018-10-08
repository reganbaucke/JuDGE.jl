module JuDGE2

push!(LOAD_PATH, ".")

using JuDGETree
using JuMP
using Gurobi


mutable struct JuDGEModel
    tree::Tree
    master::JuMP.Model
    subprob::Dict{Node,JuMP.Model}
    duals::Dict{Node,Any}
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
    jmodel.subprob::Dict{Node,JuMP.Model}()
    for n in jmodel.tree.nodes
        sp = JuMP.Model()
        jmodel.buildexpansionvariables(sp,n)
        jmodel.buildsubs(sp,n)
        jmodel.subprob[n] = sp
    end
end

function buildmaster(jmodel::JuDGEModel)
    jmodel.master = JuMP.Model()
    # jmodel.buildexpansionvariables(jmodel.master)
    # jmodel.expansioncosts(jmodel.master)
    for (key,value) in filter((key,value) -> value == :expansion,jmodel.master.ext)
        println(key)
        println(jmodel.master.objDict[key])
        println(typeof(jmodel.master.objDict[key]))
        if isa(jmodel.master.objDict[key], JuMP.Variable)
            # jmodel.master.extDict[key] = @constraint(jmodel.master, [] )
            for n in jmode.tree
                tmp = @variable(jmodel.master)
            end
            println("single")
        elseif isa(jmodel.master.objDict[key], AbstractArray)
            println("normal array")
        elseif isa(jmodel.master.objDict[key], JuMP.JuMPArray)
            println("jump array")
        end
        println("----")
    end
end

function JuDGEbuild!(jmodel::JuDGEModel)

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

macro forn(sp,x::Expr)
    tmp = "@variable(" * String(sp) * "," * String(repr(x))[3:end-1] * ", category=:Bin)"
    tmp2 = "sp.ext[:" * String(x.args[1]) * "] = :expansion"
    final = quote
        $(esc(parse(tmp)))
        $(esc(parse(tmp2)))
    end
    return final
end

macro expansion(sp,x::Expr)
    tmp = "@variable(" * String(sp) * "," * String(repr(x))[3:end-1] * ", category=:Bin)"
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

# end module
export 
JuDGEsubproblem!, JuDGEModel, solve, JuDGEexpansions!, @expansion

end