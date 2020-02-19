macro expansion(model, variable)
   ex = quote
      if !haskey($model.ext, :expansions)
         $model.ext[:expansions] = []
      end
      push!($model.ext[:expansions], @variable($model, $variable, Bin))
      $model.ext[:expansions][end]
   end
   return esc(ex)
end

macro expansionconstraint(model, name ,con)
   ex = quote
      if !haskey($model.ext, :expansionconstraints)
         $model.ext[:expansionconstraints] = []
      end
      push!($model.ext[:expansionconstraints], @constraint($model, $name ,$con))
      $model.ext[:expansions][end]
   end
   return esc(ex)
end

macro expansioncosts(model, expr)
   ex = quote
      $model.ext[:expansioncosts] = @expression($model, $expr)
   end
   return esc(ex)
end
