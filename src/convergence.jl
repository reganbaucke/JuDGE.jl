struct ConvergenceState
   obj::Float64
   ub::Float64
   lb::Float64
   abs::Float64
   rel::Float64
   time::Float64
   iter::Int64
   int::Float64
end

function Base.show(io::IO,cs::ConvergenceState)
   print(io,"")
   if cs.obj==Inf
      print(io,"          ")
   elseif cs.obj>=0
      print(io," ")
   end
   Printf.@printf(io," %e  | ",cs.obj)

   if cs.ub==Inf
      print(io,"          ")
   elseif cs.ub>=0
      print(io," ")
   end
   Printf.@printf(io,"%e ",cs.ub)

   if cs.lb==-Inf
      print(io,"         ")
   elseif cs.lb>=0
      print(io," ")
   end
   Printf.@printf(io,"%e  |  ",cs.lb)

   if cs.abs==Inf
      print(io,"          ")
   elseif cs.abs>=0
      print(io," ")
   end
   Printf.@printf(io,"%e   ",cs.abs)

   if cs.abs==Inf
      print(io,"          ")
   elseif cs.abs>=0
      print(io," ")
   end
   Printf.@printf(io,"%e  |   %e  | %9.3f  %7d",cs.rel,cs.int,cs.time,cs.iter)
end

function InitialConvergenceState()
   ConvergenceState(Inf,Inf,-Inf,Inf,Inf,0.0,0,0.0)
end

function has_converged(b::ConvergenceState, a::ConvergenceState)
   if a.abs <= b.abs || a.rel <= b.rel || a.time > b.time || a.iter > b.iter
      true
   else
      false
   end
end
