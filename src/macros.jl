macro expansion(model, variable)
   ex = quote
      if !haskey($model.ext, :expansions)
         $model.ext[:expansions] = Dict{Symbol,Any}()
      end
      if !haskey($model.ext, :forced)
         $model.ext[:forced] = Dict{Symbol,Bool}()
      end
      tmp=@variable($model, $variable, Bin)
      sym=[k for (k,v) in $model.obj_dict if v===tmp]
      $model.ext[:expansions][sym[1]]=tmp
      $model.ext[:forced][sym[1]]=false
   end
   return esc(ex)
end

macro forced_expansion(model, variable)
   ex = quote
      if !haskey($model.ext, :expansions)
         $model.ext[:expansions] = Dict{Symbol,Any}()
      end
      if !haskey($model.ext, :forced)
         $model.ext[:forced] = Dict{Symbol,Bool}()
      end
      tmp=@variable($model, $variable, Bin)
      sym=[k for (k,v) in $model.obj_dict if v===tmp]
      $model.ext[:expansions][sym[1]]=tmp
      $model.ext[:forced][sym[1]]=true
   end
   return esc(ex)
end

macro expansionconstraint(model, name ,con)
   ex = quote
      if !haskey($model.ext, :expansionconstraints)
         $model.ext[:expansionconstraints] = []
      end
      push!($model.ext[:expansionconstraints], @constraint($model, $name ,$con))
      $model.ext[:expansionconstraints][end]
   end
   return esc(ex)
end

macro expansioncosts(model, expr)
   ex = quote
      $model.ext[:expansioncosts] = @expression($model, $expr)
   end
   return esc(ex)
end

macro sp_objective(model, expr)
   ex = quote
      $model.ext[:objective]=@variable($model, obj)
      $model.ext[:objective_expr]=$expr
   end
   return esc(ex)
end
