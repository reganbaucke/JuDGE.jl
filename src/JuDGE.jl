__precompile__()

module JuDGE

using JuMP
using MathProgBase
#const JuMPVERSION = Pkg.installed("JuMP")

include("tree.jl")
include("definitions.jl")
include("deteq.jl")
include("utilities.jl")
include("output.jl")

# These interface functions just store the userinput into other functions
function JuDGEexpansions!(f!,jmodel::JuDGEModel)
    jmodel.buildexpansionvariables = f!
end

function JuDGEsubproblem!(f!,jmodel::JuDGEModel)
    jmodel.buildsubs = f!
end

function JuDGEexpansioncosts!(f!,jmodel::JuDGEModel)
    jmodel.expansioncosts = f!
end

function buildsubproblems(jmodel::JuDGEModel, s::MathProgBase.AbstractMathProgSolver)
    jmodel.subprob=Dict{Node,JuMP.Model}()
    for n in jmodel.tree.nodes
        sp = JuMP.Model(solver=s)
        jmodel.buildexpansionvariables(sp)
        jmodel.buildsubs(sp,n,sp.objDict)
        # scale objective by its probability
        sp.obj = n.p *sp.obj
        jmodel.subprob[n] = sp
    end
end

# This function updates the objective coeffecients of the expansion variables in the objective function of a particular sub problem
function updateduals(jmodel::JuDGEModel,n::Node)
    sp = jmodel.subprob[n]

    for key in keys(jmodel.mastercon[n])
        if isa(sp.objDict[key], JuMP.Variable)
            changeobjcoef!(sp.objDict[key],-getdual(jmodel.mastercon[n][key]))
        elseif isa(jmodel.mastervar[n][key], AbstractArray)
            for i in Iterators.product(indices(jmodel.mastercon[n][key])...)
                changeobjcoef!(sp.objDict[key][i],-getdual(jmodel.mastercon[n][key][i]))
            end
        elseif isa(jmodel.mastervar[n][key], JuMP.JuMPArray)
            for i in Iterators.product(jmodel.mastercon[n][key].indexsets...)
                changeobjcoef!(sp.objDict[key][i...],-getdual(jmodel.mastercon[n][key][i...]))
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
            for i in Iterators.product(indices(sp.objDict[key])...)
                obj -= getvalue(sp.objDict[key][i...])*getcoef(sp.objDict[key][i...])
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

    # grab the corresponding value of the expansion variables for use in the constraint coeffecients
    coefs = Array{Float64,1}()
    for key in keys(jmodel.mastervar[n])
        if isa(sp.objDict[key], JuMP.Variable)
            push!(coefs,getvalue(sp.objDict[key]))
        elseif isa(sp.objDict[key], AbstractArray)
            for i in Iterators.product(indices(sp.objDict[key])...)
                push!(coefs,getvalue(sp.objDict[key][i...]))
            end
        elseif isa(sp.objDict[key], JuMP.JuMPArray)
            for i in Iterators.product(sp.objDict[key].indexsets...)
                push!(coefs,getvalue(sp.objDict[key][i...]))
            end
        end
    end
    push!(coefs,1.0)

    name = "col"

    return (lb,ub,:Cont ,obj, constra, coefs, name)
end

# The main core of the JuDGE code. This code builds the master problem from the definitions of the nodal subproblems
function buildmaster(jmodel::JuDGEModel,s::MathProgBase.AbstractMathProgSolver)
    # create the jump model
    jmodel.master = JuMP.Model(solver=s)

    # initialize the dicts
    jmodel.mastervar = Dict{Node,Dict{Symbol,Any}}()
    jmodel.mastercon = Dict{Node,Dict{Symbol,Any}}()
    jmodel.covercon = Dict{Node,Any}()
    for n in jmodel.tree.nodes
        jmodel.mastervar[n] = Dict{Symbol,Any}()
        jmodel.mastercon[n] = Dict{Symbol,Any}()
    end

    # set up the variables in the master
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

                splatinto = Array{Symbol,1}()
                for i in 1:length(jmodel.mastervar[n][key].indexsets)
                    push!(splatinto,parse(alphabet[i]))
                end

                ex = Expr(:macrocall, Symbol("@constraint"), jmodel.master, parse(tmp),
                    Expr(:call,:(<=),0,Expr(:call,:sum,Expr(:generator,Expr(:ref,
                    Expr(:ref,Expr(:ref,jmodel.mastervar,:h),parse(":$key")),splatinto...),:(h = P($n))))))

                jmodel.mastercon[n][key] = eval(ex)
            end
        end
    end

    # set up the cover constraints
    for n in jmodel.tree.nodes
        jmodel.covercon[n] = @constraint(jmodel.master,0==1)
    end
end

function JuDGEbuild!(jmodel::JuDGEModel,s::MathProgBase.AbstractMathProgSolver)
    # build all the sub problems
    buildsubproblems(jmodel,s)
    for n in jmodel.tree.nodes
        jmodel.subprob[n].obj += n.p*jmodel.expansioncosts(jmodel.subprob[n],n,jmodel.subprob[n].objDict)
    end
    # build the master
    buildmaster(jmodel,s)

    # add the objective terms in the master
    for n in jmodel.tree.nodes
        jmodel.master.obj += n.p*jmodel.expansioncosts(jmodel.master,n,jmodel.mastervar[n])
    end
    jmodel.isbuilt=true

    return nothing
end

function getlowerbound(jmodel::JuDGEModel)
    lb = jmodel.master.objVal
    for n in jmodel.tree.nodes
        lb += jmodel.subprob[n].objVal - getdual(jmodel.covercon[n])
    end
    jmodel.lb = lb > jmodel.lb ? lb : jmodel.lb
    return jmodel.lb
end

function JuDGEsolve!(f,jmodel::JuDGEModel,s::MathProgBase.AbstractMathProgSolver)
    if !jmodel.isbuilt
        println("---------------------------------------")
        print("Building model...")
        JuDGEbuild!(jmodel,s)
        println("  built.")
    end

    startTime = now()
    time = Base.Dates.Millisecond(0)
    iter = 0
    ub = Inf
    lb = -Inf

    println("---------------------------------------")
    println("Solving model...")
    println("---------------------------------------")
    # perform an iteration, time is in seconds here
    keepgoing = true
    while keepgoing

        # supress the model printing stuff out
        TT = STDOUT # save original STDOUT stream
        redirect_stdout()

        (lb,ub) = iteration(jmodel)
        # for the first iteration the master problem will be infeasible
        if isnan(lb)
            lb = -Inf
        end
        if isnan(ub)
            ub = Inf
        end

        # restore io
        redirect_stdout(TT) # restore STDOUT

        iter += 1
        time = now() - startTime

        keepgoing = !f(time.value/1000,iter,lb,ub)
    end

    # just do one last solve to clean things up.
    solve(jmodel.master,suppress_warnings=true)

    println("---------------------------------------")
    println("Convergence criterion satisfied.")
    println("---------------------------------------")
end

# for use in pmap
function processnode(jmodel,n)
    updateduals(jmodel,n)
    JuMP.solve(jmodel.subprob[n])
    return buildcolumn(jmodel,n)
end

# TODO
# this function will make use of parallel processing to solve the optimisation problem
# function JuDGEpsolve!(jmodel::JuDGEModel,iter::Int64,batch::Int64)
#     if !jmodel.isbuilt
#         JuDGEbuild!(jmodel)
#     end
#     endofbatch = 0
#     for i = 1:iter
#         ind = collect(endofbatch+1:endofbatch+batch+1)
#         for i in ind
#             i = mod(i,length(jmodel.tree.nodes))
#         end
#
#         solve(jmodel.master)
#
#         columns = pmap(processnode,jmodel,jmodel.tree.nodes[ind])
#     end
# end

function iteration(jmodel::JuDGEModel)
    solved = JuMP.solve(jmodel.master,suppress_warnings=true)
    for n in jmodel.tree.nodes
        updateduals(jmodel,n)
        JuMP.solve(jmodel.subprob[n])
        # println(buildcolumn(jmodel,n))
        addcolumn(jmodel,buildcolumn(jmodel,n))
    end
    return (getlowerbound(jmodel),jmodel.master.objVal)
end

# macros for creating expansion variables
macro expansion(sp,x::Expr)
    tmp = "@variable(" * String(sp) * "," * String(repr(x))[3:end-1] * ",category=:Bin)"
    tmp2 = "sp.ext[:" * String(x.args[1]) * "] = :expansion"
    final = quote
        $(esc(parse(tmp)))
        $(esc(parse(tmp2)))
    end
    return final
end

# macros for creating expansion variables
macro expansion(sp,x::Symbol)
    tmp = "@variable(" * String(sp) * "," * String(x) * ",category=:Bin)"
    tmp2 = "sp.ext[:" * String(x) * "] = :expansion"
    final = quote
        $(esc(parse(tmp2)))
        $(esc(parse(tmp)))
    end
    return final
end


# end module
export
JuDGEsubproblem!, JuDGEModel, JuDGEsolve!, JuDGEexpansions!, @expansion, JuDGEbuild!, JuDGEexpansioncosts!, JuDGEpsolve!, JuDGEsolvedeteq!, printExpansions, printDetEqExpansions, P, Node, Tree, buildtree, getindex, getparents, getnode, getvalueDW, wholetree!, stage

end
