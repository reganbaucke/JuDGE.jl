macro expansion(model, variable)
   ex = quote
      if !haskey($model.ext, :expansions)
         $model.ext[:expansions] = Dict{Symbol,Any}()
         $model.ext[:options] = Dict{Symbol,Tuple}()
      end
      tmp=@variable($model, $variable, Bin)
      sym=[k for (k,v) in $model.obj_dict if v===tmp]
      $model.ext[:expansions][sym[1]]=tmp
      $model.ext[:options][sym[1]]=(false,0,999)
   end
   return esc(ex)
end

macro expansion(model, variable, lag)
   ex = quote
      if !haskey($model.ext, :expansions)
         $model.ext[:expansions] = Dict{Symbol,Any}()
         $model.ext[:options] = Dict{Symbol,Tuple}()
      end
      tmp=@variable($model, $variable, Bin)
      sym=[k for (k,v) in $model.obj_dict if v===tmp]
      $model.ext[:expansions][sym[1]]=tmp
      $model.ext[:options][sym[1]]=(false,$lag,999)
   end
   return esc(ex)
end

macro expansion(model, variable, lag, span)
   ex = quote
      if !haskey($model.ext, :expansions)
         $model.ext[:expansions] = Dict{Symbol,Any}()
         $model.ext[:options] = Dict{Symbol,Tuple}()
      end
      tmp=@variable($model, $variable, Bin)
      sym=[k for (k,v) in $model.obj_dict if v===tmp]
      $model.ext[:expansions][sym[1]]=tmp
      $model.ext[:options][sym[1]]=(false,$lag,$span)
   end
   return esc(ex)
end

macro shutdown(model, variable)
   ex = quote
      if !haskey($model.ext, :expansions)
         $model.ext[:expansions] = Dict{Symbol,Any}()
         $model.ext[:options] = Dict{Symbol,Tuple}()
      end
      tmp=@variable($model, $variable, Bin)
      sym=[k for (k,v) in $model.obj_dict if v===tmp]
      $model.ext[:expansions][sym[1]]=tmp
      $model.ext[:options][sym[1]]=(true,0,999)
   end
   return esc(ex)
end

macro shutdown(model, variable, lag)
   ex = quote
      if !haskey($model.ext, :expansions)
         $model.ext[:expansions] = Dict{Symbol,Any}()
         $model.ext[:options] = Dict{Symbol,Tuple}()
      end
      tmp=@variable($model, $variable, Bin)
      sym=[k for (k,v) in $model.obj_dict if v===tmp]
      $model.ext[:expansions][sym[1]]=tmp
      $model.ext[:options][sym[1]]=(true,$lag,999)
   end
   return esc(ex)
end

macro shutdown(model, variable, lag, span)
   ex = quote
      if !haskey($model.ext, :expansions)
         $model.ext[:expansions] = Dict{Symbol,Any}()
         $model.ext[:options] = Dict{Symbol,Tuple}()
      end
      tmp=@variable($model, $variable, Bin)
      sym=[k for (k,v) in $model.obj_dict if v===tmp]
      $model.ext[:expansions][sym[1]]=tmp
      $model.ext[:options][sym[1]]=(true,$lag,$span)
   end
   return esc(ex)
end

macro expansioncosts(model, expr)
   ex = quote
      $model.ext[:expansioncosts] = @expression($model, $expr)
   end
   return esc(ex)
end

macro maintenancecosts(model, expr)
   ex = quote
      $model.ext[:maintenancecosts] = @expression($model, $expr)
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
