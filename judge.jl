
using JuMP

mutable struct JuDGEModel
    tree::Tree
    master::JuMP.Model
    subprob::Array{JuMP.Model,1}
end

function capacityvariable(var)
    var.meta[:capvar] = :true
end

function capacityconstraint(con)
    con.meta[:capcon] = :true
end

function onlyoneinvestment(con)
    con.meta[:oneinvest] = :true
end



function constructmaster(m::JuMP.Model)
    # create model
    master =  Model(solver=GurobiSolver(OutputFlag=0,Method=2))

    # create capacity variables
    # capacity variables
    capvars = Symbol[]
    for key in keys(m.objDict)
        if haskey(m.objDict[key].meta,:capvar)
            @variable(master, [nodes], lowerbound = 0 ,upperbound = 1)
            push!(capvars,key)
        end
    end

    # put them in the objective function
    for key in capvars
        for var in m.objDict[key].innerArray
            if in(var,m.obj.aff.vars)
                # push!(master.obj,1)
                master.obj += 1
            end
        end
    end
    return master

end

function constructsub(m::JuMP.Model)

end
