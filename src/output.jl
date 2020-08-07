function value(jmodel::JuDGEModel, node::AbstractTree, var::Symbol, index::Int64)
   if termination_status(jmodel.master_problem) != MathOptInterface.OPTIMAL
      error("JuDGE model not solved")
   end
   return JuMP.value(jmodel.sub_problems[node][var][index])
end

function value(jmodel::JuDGEModel, node::AbstractTree, var::Symbol)
   if termination_status(jmodel.master_problem) != MathOptInterface.OPTIMAL
      error("JuDGE model not solved")
   end

   if isa(jmodel.sub_problems[node][var],JuMP.Containers.DenseAxisArray) || isa(jmodel.sub_problems[node][var],JuMP.Containers.SparseAxisArray)
       return JuMP.value.(jmodel.sub_problems[node][var])
   else
       return JuMP.value(jmodel.sub_problems[node][var])
   end
end

"""
	print_expansions(jmodel::JuDGEModel;
                    onlynonzero::Bool=true,
                    inttol=10^-9,
                    format=nothing)

Given a solved JuDGE model, this function will write the optimal capacity expansion
decisions to the REPL.

### Required Arguments
`jmodel` is the JuDGE model whose solution we wish to write to a file

### Optional Arguments
`onlynonzero` is a boolean, if set to `true` the function will only print expansions
with a non-zero value.

`inttol` is the integrality tolerance; any expansion variable value less than this
will be treated as 0, and any value greater than 1-`inttol` will be treated as 1

`format` is a function that specifies customised printing of expansion values.
See [Tutorial 2: Formatting output](@ref) for more details.
"""
function print_expansions(jmodel::JuDGEModel;onlynonzero::Bool=true,inttol=10^-9,format=nothing)
    if termination_status(jmodel.master_problem) != MathOptInterface.OPTIMAL
        error("You need to first solve the decomposed model.")
    end

    println("\nJuDGE Expansions")

    for node in collect(jmodel.tree)
        for (x,var) in jmodel.master_problem.ext[:expansions][node]
             if typeof(var) <: AbstractArray
                 if typeof(format) <: Function
                     if format_output(node,x,format(x,JuMP.value.(var)),onlynonzero,inttol)
                         continue
                     end
                 end

                 val=JuMP.value.(var)
                 for key in keys(val)
                      if !onlynonzero || val[key]>inttol
                          if typeof(val) <: Array
                              strkey=string(key)
                              strkey=replace(strkey,"CartesianIndex("=>"")
                              strkey=replace(strkey,")"=>"")
                              strkey=replace(strkey,", "=>",")
                              temp="Node "*node.name*": \""*string(x)*"["*strkey*"]\" "*string(val[key] > 1-inttol ? 1.0 : val[key])
                          else
                             temp="Node "*node.name*": \""*string(x)*"["
                             for i in 1:length(val.axes)-1
                                temp*=string(key[i])*","
                             end
                             temp*=string(key[length(val.axes)])*"]\" "*string(val[key] > 1-inttol ? 1.0 : val[key])
                         end
                         println(temp)
                     end
                 end
             else
                 if typeof(format) <: Function
                     if format_output(node,x,format(x,JuMP.value(var)),onlynonzero,inttol)
                         continue
                     end
                 end
                 if !onlynonzero || JuMP.value(var)>inttol
                     println("Node "*node.name * ": \"" * string(x) *"\" " * string(JuMP.value(var) > 1-inttol ? 1.0 : JuMP.value(var)))
                 end
             end
        end
    end
end

"""
	print_expansions(deteq::DetEqModel;
                    onlynonzero::Bool=true,
                    inttol=10^-9,
                    format=nothing)

Given a solved deterministic equivalent model, this function will write the optimal
capacity expansion decisions to the REPL.

### Required Arguments
`deteq` is the deterministic equivalent model whose solution we wish to write to a file

### Optional Arguments
`onlynonzero` is a boolean, if set to `true` the function will only print expansions
with a non-zero value.

`inttol` is the integrality tolerance; any expansion variable value less than this
will be treated as 0, and any value greater than 1-`inttol` will be treated as 1

`format` is a function that specifies customised printing of expansion values.
See [Tutorial 2: Formatting output](@ref) for more details.
"""
function print_expansions(deteq::DetEqModel;onlynonzero::Bool=true,inttol=10^-9,format=nothing)
    if termination_status(deteq.problem) != MathOptInterface.OPTIMAL
        error("You need to first solve the decomposed model.")
    end

    println("\nDeterministic Equivalent Expansions")
    for node in keys(deteq.problem.ext[:master_vars])
        for x in keys(deteq.problem.ext[:master_vars][node])
            var = deteq.problem.ext[:master_vars][node][x]
            if typeof(var)==VariableRef
                if typeof(format) <: Function
                    if format_output(node,x,format(x,JuMP.value(var)),onlynonzero,inttol)
                        continue
                    end
                end
                if !onlynonzero || JuMP.value(var)>inttol
                    name=deteq.problem.ext[:master_names][node][x]
                    println("Node "*node.name*": \""*name*"\" "*string(JuMP.value(var) > 1.0-inttol ? 1.0 : JuMP.value(var)))
                end
            elseif typeof(var) == Dict{Any,Any}
                if typeof(format) <: Function
                    temp=Dict{Any,Float64}()
                    for i in eachindex(var)
                        name=deteq.problem.ext[:master_names][node][x][i]
                        str_indices=split(name[collect(findfirst("[",name))[1]+1:collect(findfirst("]",name))[1]-1],',')
                        indices=Array{Any,1}()
                        for index in str_indices
                            if index[1]==":"
                                push!(indices,Symbol(index[2:length(index)]))
                            else
                                try
                                    push!(indices,parse(Int64,index))
                                catch
                                    push!(indices,index)
                                end
                            end
                        end
                        if length(indices)>1
                            temp[Tuple(indices)]=JuMP.value(var[i])
                        else
                            temp[indices[1]]=JuMP.value(var[i])
                        end
                    end
                    if format_output(node,x,format(x,temp),onlynonzero,inttol)
                        continue
                    end
                end

                for i in eachindex(var)
                    if !onlynonzero || JuMP.value(var[i])>inttol
                        name=deteq.problem.ext[:master_names][node][x][i]
                        println("Node "*node.name*": \""*name*"\" "*string(JuMP.value(var[i]) > 1.0-inttol ? 1.0 : JuMP.value(var[i])))
                    end
                end
            end
        end
    end
end

function format_output(node::AbstractTree,x::Symbol,exps,onlynonzero,inttol)
    if typeof(exps)==Float64 || typeof(exps)==Int64
        if !onlynonzero || abs(exps)>inttol
            println("Node "*node.name*": \""*string(x)*"\" "*string(exps))
        end
        return true
    elseif typeof(exps)==Dict{AbstractArray,Float64} || typeof(exps)==Dict{AbstractArray,Int64}
        for (key,exp) in exps
            if !onlynonzero || abs(exp)>inttol
                println("Node "*node.name*": \""*string(x)*string(key)*"\" "*string(exp))
            end
        end
        return true
    elseif typeof(exps)==Dict{Tuple,Float64} || typeof(exps)==Dict{Tuple,Int64}
        for (key,exp) in exps
            if !onlynonzero || abs(exp)>inttol
                s_key=string(key)
                s_key=replace(s_key,"("=>"")
                s_key=replace(s_key,")"=>"")
                s_key=replace(s_key,"\""=>"")
                s_key=replace(s_key,", "=>",")
                println("Node "*node.name*": \""*string(x)*"["*s_key*"]\" "*string(exp))
            end
        end
        return true
    elseif typeof(exps)==Dict{Int64,Float64} || typeof(exps)==Dict{Symbol,Float64} || typeof(exps)==Dict{String,Float64}
        for (key,exp) in exps
            if !onlynonzero || abs(exp)>inttol
                println("Node "*node.name*": \""*string(x)*"["*string(key)*"]\" "*string(exp))
            end
        end
        return true
    elseif typeof(exps)!=Nothing
        error("Formatting function must return a Float64, Dict{AbstractArray,Float64}, or nothing.")
    end
    false
end

"""
	write_solution_to_file(deteq::DetEqModel,filename::String)

Given a deterministic equivalent model and a filename, this function writes the
entire solution to a CSV.

### Required Arguments
`deteq` is the deterministic equivalent model whose solution we wish to write to a file

`filename` is the output filename
"""
function write_solution_to_file(deteq::DetEqModel,filename::String)
    if termination_status(deteq.problem) != MathOptInterface.OPTIMAL
        error("You need to first solve the decomposed model.")
    end
    file=open(filename,"w")

    println(file,"node,variable,value")

    for node in keys(deteq.problem.ext[:vars])
        for (x,var) in deteq.problem.ext[:vars][node]
            if typeof(var)==VariableRef
                println(file,string(node.name)*",\""*string(x)*"\","*string(JuMP.value(var)))
            elseif typeof(var) <: AbstractArray
                for i in eachindex(var)
                    println(file,string(node.name)*",\""*string(var[i])*"\","*string(JuMP.value(var[i])))
                end
            end
        end
    end

    for node in keys(deteq.problem.ext[:master_vars])
        for (x,var) in deteq.problem.ext[:master_vars][node]
            if typeof(var)==VariableRef
                name=deteq.problem.ext[:master_names][node][x]
                println(file,string(node.name)*",\""*string(x)*"_master\","*string(JuMP.value(var)))
            elseif typeof(var) == Dict{Any,Any}
                for i in eachindex(var)
                    name=deteq.problem.ext[:master_names][node][x][i]
                    println(file,string(node.name)*",\""*name*"_master\","*string(JuMP.value(var[i])))
                end
            end
        end
    end
    close(file)
end

"""
	write_solution_to_file(jmodel::JuDGEModel,filename::String)

Given a JuDGE model and a filename, this function writes the
entire solution to a CSV.

### Required Arguments
`jmodel` is the JuDGE model whose solution we wish to write to a file

`filename` is the output filename
"""
function write_solution_to_file(jmodel::JuDGEModel,filename::String)
    function helper(jmodel::JuDGEModel,node::AbstractTree,file::IOStream)
        vars=all_variables(jmodel.sub_problems[node])
        for v in vars
            println(file,string(node.name)*",\""*string(v)*"\","*string(JuMP.value(v)))
        end

        for (x,var) in jmodel.master_problem.ext[:expansions][node]
             if typeof(var) <: AbstractArray
                 val=JuMP.value.(var)
                 for key in keys(val)
                     temp=node.name*",\""*string(x)*"["
                     if typeof(val) <: Array
                         strkey=string(key)
                         strkey=replace(strkey,"CartesianIndex("=>"")
                         strkey=replace(strkey,")"=>"")
                         strkey=replace(strkey,", "=>",")
                         temp*=strkey
                     else
                         for i in 1:length(val.axes)-1
                            temp*=string(key[i])*","
                         end
                         temp*=string(key[length(val.axes)])
                     end
                     temp*="]_master\","*string(val[key])
                     println(file,temp)
                 end
             else
                 println(file,node.name * ",\"" * string(x) *"_master\"," * string(JuMP.value(var)))
             end
        end

        if typeof(node)==Tree
            for child in node.children
                helper(jmodel,child,file)
            end
        else
        	println(file,node.name*",\"scenario_obj\","*string(JuMP.value(jmodel.master_problem.ext[:scenprofit_var][node])))
        end
    end

    if termination_status(jmodel.master_problem) != MathOptInterface.OPTIMAL
        error("You need to first solve the decomposed model.")
    end

    file=open(filename,"w")
    println(file,"node,variable,value")
    helper(jmodel, jmodel.tree, file)

    close(file)
end
