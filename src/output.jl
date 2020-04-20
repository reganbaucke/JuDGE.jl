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

function print_expansions(jmodel::JuDGEModel;node=jmodel.tree::AbstractTree,onlynonzero::Bool=true)
    if termination_status(jmodel.master_problem) != MathOptInterface.OPTIMAL
        error("You need to first solve the decomposed model.")
    end
    # this is how you access the value of the binary expansions in the master
    for x in keys(jmodel.master_problem.ext[:expansions][node])
        var = jmodel.master_problem.ext[:expansions][node][x]
         if isa(var,Array)
             for key in keys(var)
                 if !onlynonzero || JuMP.value(var[key])>0
                     println(node.name* "_" * string(x) * "[" * string(key)* "]" * ": " * string(JuMP.value(var[key])))
                 end
             end
         elseif isa(var,JuMP.Containers.DenseAxisArray) || isa(var,JuMP.Containers.SparseAxisArray)
             val=JuMP.value.(var)
             for key in keys(val)
                  if !onlynonzero || val[key]>0
                     temp=node.name*"_"*string(x)*"["
                     for i in 1:length(val.axes)-1
                        temp*=string(key[i])*","
                     end
                     temp*=string(key[length(val.axes)])*"]:"*string(val[key])
                     println(temp)
                 end
             end
         else
             if !onlynonzero || JuMP.value(var)>0
                 println(node.name * "_" * string(x) *": " * string(JuMP.value(var)))
             end
         end
    end

    if  typeof(node)==Tree
        for child in node.children
            print_expansions(jmodel,node=child,onlynonzero=onlynonzero)
        end
    end
end

function print_expansions(deteq::DetEqModel;onlynonzero::Bool=true)
    if termination_status(deteq.problem) != MathOptInterface.OPTIMAL
        error("You need to first solve the decomposed model.")
    end
    # this is how you access the value of the binary expansions in the master

    for node in keys(deteq.problem.ext[:vars])
        for x in keys(deteq.problem.ext[:vars][node])
            if findfirst("_master",x)!=nothing
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

    println(file,"node,variable,value")

    for node in keys(deteq.problem.ext[:vars])
        for x in keys(deteq.problem.ext[:vars][node])
            if findfirst("_master",x)==nothing
                var = deteq.problem.ext[:vars][node][x][2]
                println(file,string(node.name)*",\""*x*"\","*string(JuMP.value(var)))
            end
        end
    end

    close(file)
end


function write_solution_to_file(jmodel::JuDGEModel,filename::String)
    function helper(jmodel::JuDGEModel,node::AbstractTree,file::IOStream)
        vars=all_variables(jmodel.sub_problems[node])
        for v in vars
            println(file,string(node.name)*",\""*string(v)*"\","*string(JuMP.value(v)))
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
    println(file,"node,variable,value")
    helper(jmodel, jmodel.tree, file)

    close(file)
end
