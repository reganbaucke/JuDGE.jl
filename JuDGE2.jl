module JuDGE2

push!(LOAD_PATH, ".")

using JuDGETree
using JuMP
using Gurobi

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
    covercon::Dict{Node,Any}
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
        sp = JuMP.Model(solver=GurobiSolver(OutputFlag=0))
        jmodel.buildexpansionvariables(sp)
        jmodel.buildsubs(sp,n)
        jmodel.subprob[n] = sp
    end
end

function updateduals(jmodel::JuDGEModel,n::Node)
    sp = jmodel.subprob[n]

    for key in keys(jmodel.mastercon[n])
        if isa(sp.objDict[key], JuMP.Variable)
            changeobjcoef!(sp.objDict[key],-getdual(jmodel.mastercon[n][key]))
        elseif isa(jmodel.mastercon[n], JuMP.JuMPArray)
            for i in Iterators.product(jmodel.mastercon[n].indexsets...)
                changeobjcoef!(sp.objDict[key][i],-getdual(jmodel.mastercon[n][key][i]))
            end
        end
    end
end

function addcolumn(jmodel::JuDGEModel,column)
    Variable(jmodel.master,column...)
end

function buildcolumn(jmodel::JuDGEModel,n::Node)
    sp = jmodel.subprob[n]

    lb = 0;
    ub = 1;

    # objecitve coeff is the objective value of the problem minus the terms involving the expansions
    obj = jmodel.subprob[n].objVal
    for key in keys(jmodel.mastervar[n])
        if isa(sp.objDict[key], JuMP.Variable)
            obj -= getvalue(sp.objDict[key])*getcoef(sp.objDict[key])
        elseif isa(sp.objDict[key], AbstractArray)
            for i in Iterators.product(indices(sp.objDict[key]))
                obj -= getvalue(sp.objDict[key][i])*getcoef(sp.objDict[key][i])
            end
        elseif isa(sp.objDict[key], JuMP.JuMPArray)
            for i in Iterators.product(sp.objDict[key].indexsets...)
                obj -= getvalue(sp.objDict[key][i...])*getcoef(sp.objDict[key][i...])
            end
        end
    end

    # constraints collection is are the expansion variables at this node
    # println(jmodel.mastervar)
    constra = []
    for key in keys(jmodel.mastervar[n])
        if isa(jmodel.mastervar[n][key], JuMP.Variable)
            push!(constra,jmodel.mastercon[n][key])
        elseif isa(jmodel.mastercon[n][key], JuMP.JuMPArray)
            for i in Iterators.product(jmodel.mastercon[n][key].indexsets...)
                push!(constra,jmodel.mastercon[n][key][i...])
            end
        elseif isa(jmodel.mastercon[n][key], AbstractArray)
            for i in Iterators.product(indices(jmodel.mastercon[n][key])...)
                push!(constra,jmodel.mastercon[n][key][i...])
            end
        end
    end
    push!(constra,jmodel.covercon[n])

    coefs = Array{Float64,1}()
    for key in keys(jmodel.mastervar[n])
        if isa(sp.objDict[key], JuMP.Variable)
            push!(coefs,getvalue(sp.objDict[key]))
        elseif isa(sp.objDict[key], AbstractArray)
            for i in Iterators.product(indices(sp.objDict[key]))
                push!(coefs,getvalue(sp.objDict[key][i]))
            end
        elseif isa(sp.objDict[key], JuMP.JuMPArray)
            for i in Iterators.product(sp.objDict[key].indexsets...)
                push!(coefs,getvalue(sp.objDict[key][i...]))
            end
        end
    end
    push!(coefs,1.0)

    # contr = [ jmodel.master.objDict[:pi][n]; hello.master.objDict[:mu][n]]
    # colcoef = [getvalue(sp[:z]),1]
    name = "James"

    return (lb,ub,:Cont ,obj, constra, coefs, name)
end

function changeobjcoef!(var::JuMP.Variable,coef::Float64)
    # find model
    m = var.m
    pos = findfirst(m.obj.aff.vars,var)
    if pos == 0
        m.obj += var*coef
    else
        m.obj.aff.coeffs[pos] = coef
    end
end

function getcoef(var::JuMP.Variable)
    # find model
    m = var.m
    pos = findfirst(m.obj.aff.vars,var)
    if pos == 0
        return 0.0
    else
        return m.obj.aff.coeffs[pos]
    end
end

function buildmaster(jmodel::JuDGEModel)
    # create the jump model
    jmodel.master = JuMP.Model(solver=GurobiSolver(OutputFlag=0))

    # initialize the dicts
    jmodel.mastervar = Dict{Node,Dict{Symbol,Any}}()
    jmodel.mastercon = Dict{Node,Dict{Symbol,Any}}()
    jmodel.covercon = Dict{Node,Any}()
    for n in jmodel.tree.nodes
        jmodel.mastervar[n] = Dict{Symbol,Any}()
        jmodel.mastercon[n] = Dict{Symbol,Any}()
    end

    # set up the variables
    sp = jmodel.subprob[jmodel.tree.root]
    for (key,value) in filter((key,value) -> value == :expansion,sp.ext)
        if isa(sp.objDict[key], JuMP.Variable)
            for n in jmodel.tree.nodes
                jmodel.mastervar[n][key] = @variable(jmodel.master,lowerbound = 0, upperbound = 1)
            end
        elseif isa(sp.objDict[key], AbstractArray)
            for n in jmodel.tree.nodes
                # have to do it this hacker way because the only nice way to make variables is with @variable
                ex = Expr(:macrocall, Symbol("@variable"), jmodel.master, Expr(:vect, indices(sp.objDict[key])...),:(lowerbound = 0), :(upperbound = 1) )
                jmodel.mastervar[n][key] = eval(ex)
            end
        elseif isa(sp.objDict[key], JuMP.JuMPArray)
            for n in jmodel.tree.nodes
                # have to do it this hacker way because the only nice way to make variables is with @variable
                ex = Expr(:macrocall, Symbol("@variable"), jmodel.master, Expr(:vect, sp.objDict[key].indexsets...), :(lowerbound = 0), :(upperbound = 1))
                jmodel.mastervar[n][key] = eval(ex)
            end
        end
    end

    # set up the constraints
    for (key,value) in filter((key,value) -> value == :expansion,sp.ext)
        if isa(sp.objDict[key], JuMP.Variable)
            for n in jmodel.tree.nodes
                jmodel.mastercon[n][key] = @constraint(jmodel.master, 0 <= sum(jmodel.mastervar[h][key] for h in P(n)))
            end
        elseif isa(sp.objDict[key], AbstractArray) || isa(sp.objDict[key],JuMP.JuMPArray)
            
            for n in jmodel.tree.nodes
                # we need this to define the counters which loop of the index sets
                tmp = collect('a':'z')
                alphabet = Array{String,1}()
                for i in 1:length(tmp)
                    push!(alphabet,string(tmp[i]))
                end

                tmp = "["
                for (i,set) in enumerate(jmodel.mastervar[n][key].indexsets)
                    tmp *=  "$(alphabet[i])" *" in "  * repr(set)
                    if i != length(jmodel.mastervar[n][key].indexsets)
                        tmp *= ","
                    end
                end
                tmp *= "]"

                tmp2 = "0 <= sum(jmodel.mastervar[h][key]["
                for (i,set) in enumerate(jmodel.mastervar[n][key].indexsets)
                    tmp2 *=  "a"
                    if i != length(jmodel.mastervar[n][key].indexsets)
                        tmp2 *= ","
                    end
                end
                tmp2 *= "] for h in P(n))"

                # println(key)
                # println(Symbol(key))
                # something = Expr(:ref, jmodel.mastervar[jmodel.tree.root], parse(":$key") )
                # println(something)
                # eval(something)

                splatinto = Array{Symbol,1}()
                for i in 1:length(jmodel.mastervar[n][key].indexsets)
                    push!(splatinto,parse(alphabet[i]))
                end

                ex = Expr(:macrocall, Symbol("@constraint"), jmodel.master, parse(tmp),
                    Expr(:call,:(<=),0,Expr(:call,:sum,Expr(:generator,Expr(:ref,
                    Expr(:ref,Expr(:ref,jmodel.mastervar,:h),parse(":$key")),splatinto...),:(h = P($n))))))

                jmodel.mastercon[n][key] = eval(ex)
            end
        elseif isa(sp.objDict[key], JuMP.JuMPArray)
            for n in jmodel.tree.nodes
                # have to do it this hacker way because the only nice way to make variables is with @variable
                # ex = Expr(:macrocall, Symbol("@variable"), sp, Expr(:vect, sp.objDict[key].indexsets...))
                # jmodel.mastervar[n][key] = eval(ex)
            end
        end
    end
    # set up the cover constraints
    for n in jmodel.tree.nodes
        jmodel.covercon[n] = @constraint(jmodel.master,0==1)
    end

end

function JuDGEbuild!(jmodel::JuDGEModel)
    # build all the sub problems
    buildsubproblems(jmodel)

    # build the master
    buildmaster(jmodel)

    # add the objective terms in the master
    for n in jmodel.tree.nodes
        jmodel.master.obj += jmodel.expansioncosts(jmodel.master,jmodel.mastervar[n],n)
    end
end


function JuDGEsolve!(jmodel::JuDGEModel,iter::Int64)
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
            updateduals(jmodel,n)
        end
        JuMP.solve(jmodel.subprob[n])
        addcolumn(jmodel,buildcolumn(jmodel,n))
    end
end

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
JuDGEsubproblem!, JuDGEModel, JuDGEsolve!, JuDGEexpansions!, @expansion, JuDGEbuild!, JuDGEexpansioncosts!

end