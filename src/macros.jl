macro judge_var(model, variable, class, aargs, aakws)
   lag=0
   span=1000
   initial=0.0
   ub=nothing
   lb=nothing
   penalty=nothing
   state_name=:nothing

   vartype=:Con

   if length(aargs)==0
      push!(aargs,:Con)
   end

   if length(aargs)!=1
      ex = quote
         error("@"*string($class)*" macro takes at most three positional arguments")
      end
      return ex
   elseif aargs[1] ∉ [:Con,:Bin,:Int]
      ex = quote
         error("Optional third positional argument must be \'Bin\' or \'Int\'")
      end
      return ex
   end

   for (a,b) in aakws
      if a==:lag
         lag=b
      elseif a==:duration
         span=b
      elseif a==:initial
         initial=b
      elseif a==:lb
         lb=b
      elseif a==:ub
         ub=b
      elseif a==:penalty
         penalty=b
      elseif a==:state_name
         state_name=b
      else
         ex = quote
            error("Invalid keyword argument for @"*string($class)*" macro")
         end
         return ex
      end
   end

   if class==:(:state)
      if lag!=0
         @warn("'lag' keyword has been ignored")
         lag=0
      end
      if span!=1000
         @warn("'duration' keyword has been ignored")
      end
      span=1
   elseif state_name!=:nothing
      @warn("'state_name' keyword has been ignored")
      state_name=:nothing
   end

   if !(class==:state || class==:enforced) && penalty!=nothing
      @warn("'penalty' keyword has been ignored")
      penalty=nothing
   end

   tmp=nothing

   ex = quote
      if !haskey($model.ext, :expansions)
         $model.ext[:expansions] = Dict{Symbol,Any}()
         $model.ext[:options] = Dict{Symbol,Tuple}()
      end
      if length($aargs)==1
         if $aargs[1]==:Con
            tmp=@variable($model, $variable)
         elseif $aargs[1]==:Bin
            tmp=@variable($model, $variable, Bin)
         else
            tmp=@variable($model, $variable, Int)
         end
      end

      if $(Meta.quot(state_name))==:nothing
         sym=[k for (k,v) in $model.obj_dict if v===tmp][1]
         $model.ext[:expansions][sym]=tmp
         $model.ext[:options][sym]=($class,$lag,$span,$aargs[1],$lb,$ub,$initial,$penalty)
      else
         $model.ext[:expansions][$(Meta.quot(state_name))]=tmp
         $model.ext[:options][$(Meta.quot(state_name))]=($class,$lag,$span,$aargs[1],$lb,$ub,$initial,$penalty)
         $state_name=tmp
      end

   end
   return esc(ex)
end

"""
	expansion(model, variable, vartype, lag=0, span=1000)

Defines an expansion variable `variable` within a subproblem `model`. Note that all subproblems must have the same set of expansion variables.

### Required Arguments
`model` is the JuDGE subproblem that we are adding the expansion variable to

`variable` is the name of the variable being created, this will be automatically set to be binary; follows JuMP syntax if defining a set of variables.

### Optional Arguments
`lag` is the number of nodes in the scenario between an expansion being decided, and it becoming available

`span` is the number of consecutive nodes in the scenario over which an expansion is available

### Examples
    @expansion(model, expand[1:5]) #defines an array of 5 variables with no lag, and unlimited lifespan
    @expansion(model, expand[1:5,1:2], 1) #defines a matrix of 10 variables with a lag of 1, and unlimited lifespan
    @expansion(model, expand, 0, 2) #defines a single variable with a lag of 0, and a lifespan of 2
"""
macro expansion(model, variable, args...)
   aargs = []
   aakws = Pair{Symbol,Any}[]
   for el in args
      if Meta.isexpr(el, :(=))
         push!(aakws, Pair(el.args...))
      else
         push!(aargs, el)
      end
   end
   ex = quote
      JuDGE.@judge_var($model, $variable, :expansion, $aargs, $aakws)
   end

   return esc(ex)
end

"""
	shutdown(model, variable, vartype, lag=0, span=1000)

Defines an shutdown variable `variable` within a subproblem `model`. Note that all subproblems must have the same set of shutdown variables.

### Required Arguments
`model` is the JuDGE subproblem that we are adding the shutdown variable to

`variable` is the name of the variable being created, this will be automatically set to be binary; follows JuMP syntax if defining a set of variables.

### Optional Arguments
`lag` is the number of nodes in the scenario between an shutdown being decided, and it becoming unavailable

`span` is the number of consecutive nodes in the scenario over which the shutdown will last

### Examples
    @shutdown(model, shut[1:5]) #defines an array of 5 variables with no lag, and unlimited duration
    @shutdown(model, shut[1:5,1:2], 1) #defines a matrix of 10 variables with a lag of 1, and unlimited duration
    @shutdown(model, shut, 0, 2) #defines a single variable with a lag of 0, and a lifespan of 2
"""
macro shutdown(model, variable, args...)
   aargs = []
   aakws = Pair{Symbol,Any}[]
   for el in args
      if Meta.isexpr(el, :(=))
         push!(aakws, Pair(el.args...))
      else
         push!(aargs, el)
      end
   end
   ex = quote
      JuDGE.@judge_var($model, $variable, :shutdown, $aargs, $aakws)
   end

   return esc(ex)
end

macro enforced(model, variable, args...)
   aargs = []
   aakws = Pair{Symbol,Any}[]
   for el in args
      if Meta.isexpr(el, :(=))
         push!(aakws, Pair(el.args...))
      else
         push!(aargs, el)
      end
   end
   ex = quote
      JuDGE.@judge_var($model, $variable, :enforced, $aargs, $aakws)
   end

   return esc(ex)
end

macro state(model, variable, args...)
   aargs = []
   aakws = Pair{Symbol,Any}[]
   for el in args
      if Meta.isexpr(el, :(=))
         push!(aakws, Pair(el.args...))
      else
         push!(aargs, el)
      end
   end
   ex = quote
      JuDGE.@judge_var($model, $variable, :state, $aargs, $aakws)
   end

   return esc(ex)
end

"""
	capitalcosts(model, expr)

Defines a linear expression specifying the capital cost of expansions and shutdowns at the current node

### Required Arguments
`model` is the JuDGE subproblem corresponding to the node in the scenario tree that we are adding specifying the costs for

`expr` is an `AffExpr` which gives the total cost of choosing expansion and shutdown variables at the current node

### Example
    @capitalcosts(model, sum(expand[i]*cost[node][i] for i in 1:5))
"""
macro capitalcosts(model, expr)
   ex = quote
      #$model.ext[:capitalcosts] = @expression($model, $expr)
      $model.ext[:capitalcosts] = Dict()
      if typeof($expr)==AffExpr
         for (term,coef) in $expr.terms
            $model.ext[:capitalcosts][term]=coef
         end
         $model.ext[:capitalcosts][:constant]=$expr.constant
      elseif typeof($expr)==VariableRef
         $model.ext[:capitalcosts][$expr]=1.0
         $model.ext[:capitalcosts][:constant]=0.0
      elseif typeof($expr)==Float64
         $model.ext[:capitalcosts][:constant]=$expr
      end
   end
   return esc(ex)
end

"""
	ongoingcosts(model, expr)

Defines a linear expression specifying the ongoing costs of expansions and shutdowns available at the current node

### Required Arguments
`model` is the JuDGE subproblem corresponding to the node in the scenario tree that we are adding specifying the costs for

`expr` is an `AffExpr` which gives the ongoing cost of expansions and shutdowns available at the current node

### Example
    @ongoingcosts(model, sum(expand[i]*ongoingcosts[node][i] for i in 1:5))
"""
macro ongoingcosts(model, expr)
   ex = quote
#      $model.ext[:ongoingcosts] = @expression($model, $expr)
      $model.ext[:ongoingcosts] = Dict()
      if typeof($expr)==AffExpr
         for (term,coef) in $expr.terms
            $model.ext[:ongoingcosts][term]=coef
         end
         $model.ext[:ongoingcosts][:constant]=$expr.constant
      elseif typeof($expr)==VariableRef
         $model.ext[:ongoingcosts][$expr]=1.0
         $model.ext[:ongoingcosts][:constant]=0.0
      elseif typeof($expr)==Float64
         $model.ext[:ongoingcosts][:constant]=$expr
      end
   end
   return esc(ex)
end
