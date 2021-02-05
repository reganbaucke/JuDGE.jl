macro expansion(model, variable)
   ex = quote
      if !haskey($model.ext, :expansions)
         $model.ext[:expansions] = Dict{Symbol,Any}()
         $model.ext[:options] = Dict{Symbol,Tuple}()
      end
      tmp=@variable($model, $variable, Bin)
      sym=[k for (k,v) in $model.obj_dict if v===tmp]
      $model.ext[:expansions][sym[1]]=tmp
      $model.ext[:options][sym[1]]=(false,0,999,false)
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
      $model.ext[:options][sym[1]]=(false,$lag,999,false)
   end
   return esc(ex)
end

"""
	expansion(model, variable, lag, span)

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
macro expansion(model, variable, lag, span)
   ex = quote
      if !haskey($model.ext, :expansions)
         $model.ext[:expansions] = Dict{Symbol,Any}()
         $model.ext[:options] = Dict{Symbol,Tuple}()
      end
      tmp=@variable($model, $variable, Bin)
      sym=[k for (k,v) in $model.obj_dict if v===tmp]
      $model.ext[:expansions][sym[1]]=tmp
      $model.ext[:options][sym[1]]=(false,$lag,$span,false)
   end
   return esc(ex)
end

macro expansion_cont(model, variable, lag, span)
   ex = quote
      if !haskey($model.ext, :expansions)
         $model.ext[:expansions] = Dict{Symbol,Any}()
         $model.ext[:options] = Dict{Symbol,Tuple}()
      end
      tmp=@variable($model, $variable)
      sym=[k for (k,v) in $model.obj_dict if v===tmp]
      $model.ext[:expansions][sym[1]]=tmp
      $model.ext[:options][sym[1]]=(false,$lag,$span,true)
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
      $model.ext[:options][sym[1]]=(true,0,999,false)
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
      $model.ext[:options][sym[1]]=(true,$lag,999,false)
   end
   return esc(ex)
end

"""
	shutdown(model, variable, lag, span)

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
      for (term,coef) in $expr.terms
         $model.ext[:capitalcosts][term]=coef
      end
      $model.ext[:capitalcosts][:constant]=$expr.constant
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
      for (term,coef) in $expr.terms
         $model.ext[:ongoingcosts][term]=coef
      end
      $model.ext[:ongoingcosts][:constant]=$expr.constant
   end
   return esc(ex)
end
