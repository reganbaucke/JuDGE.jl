# Prints the expansions occurring at each node from the decomposed model
function print_expansions(jmodel::JuDGEModel;node=jmodel.tree.root::Node,onlynonzero::Bool=true)
    if !jmodel.isbuilt
        error("You need to first solve the decomposed model.")
    end
    # this is how you access the value of the binary expansions in the master
    for x in keys(jmodel.subprob[node].ext)
        var = jmodel.mastervar[node][x]
        # print value of variable (if non-zero)
        if isa(var,JuMP.JuMPArray)
            for key in keys(var)
                if !onlynonzero || getvalue(var[key...])>0
                    println(string(node) * "_" * string(x) * "[" * string(key)[2:length(string(key))-1] * "]" * ": " * string(getvalue(var[key...])))
                end
            end
        else
            println(string(node) * "_" * string(x) *": " * string(getvalue(var)))
        end
    end

    if !(length(node.children) == 0)
        for i in 1:length(node.children)
            print_expansions(jmodel,node=node.children[i],onlynonzero=onlynonzero)
        end
    end
end

# Prints the expansions occurring at each node from the deterministic equivalent
# model.
function printDetEqExpansions(jmodel::JuDGEModel;onlynonzero::Bool=true)
    if !jmodel.isbuiltdeteq
        error("You need to first solve the deterministic equivalent.")
    end
    m=jmodel.deteq

    for x in keys(jmodel.subprob[jmodel.tree.root].ext)
        expansion=string(x)
        for index in 1:m.numCols
            if contains(getname(Variable(m,index)),expansion)
                if getvalue(Variable(m,index)) > 0
                    print(getname(Variable(m,index)) * ": " * string(getvalue(Variable(m,index))) * '\n')
                end
            end
        end
    end
end
