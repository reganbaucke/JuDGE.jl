macro expansion(model, variable)
   ex = quote
      if !haskey($model.ext, :expansions)
         $model.ext[:expansions] = Dict{Any,Any}()
      end

      tmp=@variable($model, $variable, Bin)
      sym=[k for (k,v) in $model.obj_dict if v===tmp]
      $model.ext[:expansions][sym[1]]=tmp
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
