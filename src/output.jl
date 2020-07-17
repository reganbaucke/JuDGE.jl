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

function print_expansions(jmodel::JuDGEModel;onlynonzero::Bool=true,inttol=10^-9)
    if termination_status(jmodel.master_problem) != MathOptInterface.OPTIMAL
        error("You need to first solve the decomposed model.")
    end

    println("\nJuDGE Expansions")

    for node in collect(jmodel.tree)
        for (x,var) in jmodel.master_problem.ext[:expansions][node]
             if typeof(var) <: AbstractArray
                 val=JuMP.value.(var)
                 for key in keys(val)
                      if !onlynonzero || val[key]>inttol
                          if typeof(val) <: Array
                              strkey=string(key)
                              strkey=replace(strkey,"CartesianIndex("=>"")
                              strkey=replace(strkey,")"=>"")
                              strkey=replace(strkey,", "=>",")
                              temp=node.name*"_"*string(x)*"["*strkey*"]:"*string(val[key])
                          else
                             temp=node.name*"_"*string(x)*"["
                             for i in 1:length(val.axes)-1
                                temp*=string(key[i])*","
                             end
                             temp*=string(key[length(val.axes)])*"]:"*string(val[key])
                         end
                         println(temp)
                     end
                 end
             else
                 if !onlynonzero || JuMP.value(var)>0
                     println(node.name * "_" * string(x) *": " * string(JuMP.value(var)))
                 end
             end
        end
    end
end

function print_expansions(deteq::DetEqModel;onlynonzero::Bool=true)
    if termination_status(deteq.problem) != MathOptInterface.OPTIMAL
        error("You need to first solve the decomposed model.")
    end

    println("\nDeterministic Equivalent Expansions")
    for node in keys(deteq.problem.ext[:vars])
        for x in keys(deteq.problem.ext[:vars][node])
            if typeof(x)==String && findfirst("_master",x)!=nothing
                var = deteq.problem.ext[:vars][node][x]
                if !onlynonzero || JuMP.value(var)>0
                    println(node.name * "_" * replace(x,"_master"=>"") *": " * string(JuMP.value(var)))
                end
            end
        end
    end
end

function write_solution_to_file(deteq::DetEqModel,filename::String)
    if termination_status(deteq.problem) != MathOptInterface.OPTIMAL
        error("You need to first solve the decomposed model.")
    end
    file=open(filename,"w")

    println(file,"node,variable,value,obj_coeff")

    for node in keys(deteq.problem.ext[:vars])
        for (x,var) in deteq.problem.ext[:vars][node]
            println(file,string(node.name)*",\""*string(x)*"\","*string(JuMP.value(var))*","*string(objcoef(var)))
        end
    end

    close(file)
end

function write_solution_to_file(jmodel::JuDGEModel,filename::String)
    function helper(jmodel::JuDGEModel,node::AbstractTree,file::IOStream)
        vars=all_variables(jmodel.sub_problems[node])
        for v in vars
            println(file,string(node.name)*",\""*string(v)*"\","*string(JuMP.value(v))*","*string(objcoef(v)))
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
                     temp*="]_master\","*string(val[key])*","*string(objcoef(var[key]))
                     println(file,temp)
                 end
             else
                 println(file,node.name * ",\"" * string(x) *"_master\"," * string(JuMP.value(var))*","*string(objcoef(var)))
             end
        end

        if typeof(node)==Tree
            for child in node.children
                helper(jmodel,child,file)
            end
        end
    end

    if termination_status(jmodel.master_problem) != MathOptInterface.OPTIMAL
        error("You need to first solve the decomposed model.")
    end

    file=open(filename,"w")
    println(file,"node,variable,value,obj_coeff")
    helper(jmodel, jmodel.tree, file)

    close(file)
end
