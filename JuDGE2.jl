module JuDGE2

push!(LOAD_PATH, ".")

using JuDGETree
using JuMP
# using Gurobi

function P(n::Node)
    list = Node[]
    list = getparents(n)
    push!(list,n)
end

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
    jmodel.mastervar = Dict{Node,Dict{Symbol,Any}}()
    jmodel.mastercon = Dict{Node,Dict{Symbol,Any}}()
    for n in jmodel.tree.nodes
        jmodel.mastervar[n] = Dict{Symbol,Any}()
        jmodel.mastercon[n] = Dict{Symbol,Any}()
    end

    # set up the variables
    sp = jmodel.subprob[jmodel.tree.root]
    for (key,value) in filter((key,value) -> value == :expansion,sp.ext)
        if isa(sp.objDict[key], JuMP.Variable)
            for n in jmodel.tree.nodes
                jmodel.mastervar[n][key] = @variable(jmodel.master)
            end
        elseif isa(sp.objDict[key], AbstractArray)
            for n in jmodel.tree.nodes
                # have to do it this hacker way because the only nice way to make variables is with @variable
                ex = Expr(:macrocall, Symbol("@variable"), jmodel.master, Expr(:vect, indices(sp.objDict[key])...) )
                jmodel.mastervar[n][key] = eval(ex)
            end
        elseif isa(sp.objDict[key], JuMP.JuMPArray)
            for n in jmodel.tree.nodes
                # have to do it this hacker way because the only nice way to make variables is with @variable
                ex = Expr(:macrocall, Symbol("@variable"), jmodel.master, Expr(:vect, sp.objDict[key].indexsets...))
                jmodel.mastervar[n][key] = eval(ex)
            end
        end
    end

    # set up the constraints
    for (key,value) in filter((key,value) -> value == :expansion,sp.ext)
        println(key)
        if isa(sp.objDict[key], JuMP.Variable)
            for n in jmodel.tree.nodes
                jmodel.mastercon[n][key] = @constraint(jmodel.master, 0 <= sum(jmodel.mastervar[h][key] for h in P(n)))
            end
        elseif isa(sp.objDict[key], AbstractArray)
            for n in jmodel.tree.nodes
                # if n == jmodel.tree.root
                    tmp = "@constraint(jmodel.master,["
                    for (i,set) in enumerate(jmodel.mastervar[n][key].indexsets)
                        tmp *=  string(i) *" in "  * repr(set)
                        if i != length(jmodel.mastervar[n][key].indexsets)
                            tmp *= ","
                        end
                    end
                    tmp *= "], 0 <= sum( jmodel.mastervar[h][key]["

                    for (i,set) in enumerate(jmodel.mastervar[n][key].indexsets)
                        tmp *=  string(i)
                        if i != length(jmodel.mastervar[n][key].indexsets)
                            tmp *= ","
                        end
                    end
                    tmp *= "] for h in P(n)))"

                    println(tmp)
                    eval(parse(tmp))
                # end
            end
        elseif isa(sp.objDict[key], JuMP.JuMPArray)
            for n in jmodel.tree.nodes
                # have to do it this hacker way because the only nice way to make variables is with @variable
                # ex = Expr(:macrocall, Symbol("@variable"), sp, Expr(:vect, sp.objDict[key].indexsets...))
                # jmodel.mastervar[n][key] = eval(ex)
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


# end module
export 
JuDGEsubproblem!, JuDGEModel, solve, JuDGEexpansions!, @expansion, JuDGEbuild!

end