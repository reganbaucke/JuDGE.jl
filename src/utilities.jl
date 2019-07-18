# helper function which changes the objective function coeffecient of a particular variable
# this really should be a jump standard
function changeobjcoef!(var::JuMP.Variable,coef::Float64)
    # find model
    m = var.m
    pos = findfirst(m.obj.aff.vars,var)
    if pos == 0
        m.obj += var*coef
    else
        m.obj.aff.coeffs[pos] = coef

        while (pos=findnext(m.obj.aff.vars,var,pos+1))!=0
            m.obj.aff.coeffs[pos]=0
        end
    end
end

# helper function which fetches the current objective value coef for a given variable
function getcoef(var::JuMP.Variable)
    # find model
    m = var.m
    pos = findfirst(m.obj.aff.vars,var)
    if pos == 0
        return 0.0
    else
        coeff=m.obj.aff.coeffs[pos]
        while (pos=findnext(m.obj.aff.vars,var,pos+1))!=0
            coeff+=m.obj.aff.coeffs[pos]
        end
        return coeff
    end
end

# get values of a variable from a node
function getvalueDW(jmodel::JuDGEModel, node::Node, var::Symbol)
    if jmodel.isbuilt
        return getvalue(jmodel.subprob[node][var])
    else
        println("Decomposition model not built.")
    end
end

function getvalueDW(jmodel::JuDGEModel, indices::Array{Int64,1}, var::Symbol)
    if jmodel.isbuilt
        return getvalue(jmodel.subprob[getnode(jmodel.tree,indices)][var])
    else
        println("Decomposition model not built.")
    end
end

function JuDGEwriteLP(jmodel::JuDGEModel, LP::Symbol, filename::String)
    if LP==:DetEq
        if jmodel.isbuiltdeteq
            f = open(filename, "w")
            print(f, jmodel.deteq)
            close(f)
        else
            println("Deterministic equivalent model has not yet been built")
        end
    elseif LP==:Master
        if jmodel.isbuilt
            f = open(filename, "w")
            print(f, jmodel.master)
            close(f)
        else
            println("Danzig-wolfe model has not yet been built")
        end
    elseif LP==:Subproblems
        if jmodel.isbuilt
            queue=[jmodel.tree.root]
            while size(queue)[1]>0
                node=pop!(queue);
                f = open(filename * "_" * replace(SubString(string(node),6,length(string(node))-1),", ","_") * ".lp", "w")
                print(f, jmodel.subprob[node])
                close(f)
                for i in node.children
                    push!(queue,i)
                end
            end
        else
            println("Danzig-wolfe model has not yet been built")
        end
    else
        println("Accepted symbols for LP are: :DetEq, :Master, :Subproblems")
    end
end
