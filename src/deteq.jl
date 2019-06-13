# Builds a deterministic equivalent by joining the subproblems into
# a single MIP, then adding the constraints that there can only be
# a single investment along any path in the tree.
function JuDGEbuilddeteq(jmodel::JuDGEModel, s::MathProgBase.AbstractMathProgSolver)
    jmodel.deteq = JuMP.Model(solver=s)

    # Always rebuild the subproblems so they have the correct
    # costs. Can be easily made more efficient.
    JuDGEbuild!(jmodel,s)

    offset=0
    ex=Expr(:macrocall, Symbol("@objective"),jmodel.deteq, :(Min), 0)
    eval(ex)
    history=Dict{Symbol,Array{Int64,1}}()

    queue=[]

    push!(queue,jmodel.tree.root)

    #for n in jmodel.tree.nodes
    while length(queue)>0
        n=pop!(queue)
        println("Processing subproblem " * string(n))
        for c in n.children
            push!(queue,c)
        end
        expansions=[]
        sp = jmodel.subprob[n]

        for i in 1:sp.numCols
          v=Variable(sp,i)
          lastchr=findin(string(v),"[")[1]
          if length(lastchr)==0
              lastchr=length(string(v))+1
          end
          #if sp.ext[Symbol(string(v)[1:lastchr-1])] == :expansion
          if Symbol(string(v)[1:lastchr-1]) in keys(sp.ext)
              if !(Symbol(v) in keys(history))
                  history[Symbol(v)]=[]
              end
              push!(history[Symbol(v)],offset+i)
              push!(expansions,i)
          end
          name=getname(v)
          ex=Expr(:macrocall, Symbol("@variable"),jmodel.deteq,Symbol(replace(string(n),r"[|]"=>"!") * "_" * name))
          eval(ex)
          jmodel.deteq.colCat[jmodel.deteq.numCols]=sp.colCat[i]
          jmodel.deteq.colLower[jmodel.deteq.numCols]=sp.colLower[i]
          jmodel.deteq.colUpper[jmodel.deteq.numCols]=sp.colUpper[i]

          changeobjcoef!(Variable(jmodel.deteq,jmodel.deteq.numCols),getcoef(v))
        end

        if length(n.children)==0
            for j in 1:numCols
                vars=[]
                coeffs=[]
                v=Variable(jmodel.deteq,offset+j)
                v2=Variable(sp,j)
                if j in expansions
                    for k in history[Symbol(v2)]
                        v3=Variable(jmodel.deteq,k)
                        tmpstr=string(v3)[1:findin(string(v3),"_")[1]-2]
                        if contains(string(v),tmpstr)
                            push!(vars,v3)
                            push!(coeffs,1)
                        end
                    end
                    LHS = AffExpr(vars, coeffs, 0)
                    @constraint(jmodel.deteq,LHS<=1)
                end
            end
        end

        concoeffs=JuMP.prepConstrMatrix(sp)
        RHSs=JuMP.prepConstrBounds(sp)
        numRows,numCols=size(concoeffs)

        for i in 1:numRows
            vars=[]
            coeffs=[]
            for j in 1:numCols
                v=Variable(jmodel.deteq,offset+j)
                v2=Variable(sp,j)
                if concoeffs[i,j]!=0
                    if j in expansions
                        for k in history[Symbol(v2)]
                            v3=Variable(jmodel.deteq,k)
                            tmpstr=string(v3)[1:findin(string(v3),"_")[1]-2]
                            if contains(string(v),tmpstr)
                                push!(vars,v3)
                                push!(coeffs,concoeffs[i,j])
                            end
                        end
                    else
                        push!(vars,v)
                        push!(coeffs,concoeffs[i,j])
                    end
                end
            end
            LHS = AffExpr(vars, coeffs, 0)
            if RHSs[1][i]==RHSs[2][i]
                @constraint(jmodel.deteq,LHS==RHSs[1][i])
            elseif RHSs[1][i]==-Inf
                @constraint(jmodel.deteq,LHS<=RHSs[2][i])
            elseif RHSs[2][i]==Inf
                @constraint(jmodel.deteq,LHS>=RHSs[1][i])
            else
                @constraint(jmodel.deteq,RHSs[1][i]<=LHS<=RHSs[2][i])
            end
        end
        offset+=numCols
    end

    jmodel.isbuiltdeteq=true

    return nothing
end

# Function called by the user to solve the deterministic equivalent
function JuDGEsolvedeteq!(jmodel::JuDGEModel,s::MathProgBase.AbstractMathProgSolver)
    if !jmodel.isbuiltdeteq
        println("---------------------------------------")
        print("Building deterministic equivalent formulation...")
        JuDGEbuilddeteq(jmodel,s)
        println("  built.")
        jmodel.isbuiltdeteq=true
    end

    status = solve(jmodel.deteq)
    if status == :Optimal
        println("  solved.")
    else
        println("  not solved: " * string(status))
    end
end
