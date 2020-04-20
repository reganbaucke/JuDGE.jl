struct ConvergenceState
   ub::Float64
   lb::Float64
   abs::Float64
   rel::Float64
   time::Float64
   iter::Int64
   int::Float64
end

function Base.show(io::IO,cs::ConvergenceState)
   temp=""

   if cs.ub==Inf
      print(io,"           ")
   elseif cs.ub<0
      print(io," ")
   end
   Printf.@printf(io,"%e ",cs.ub)

   if cs.lb==-Inf
      print(io,"         ")
   elseif cs.ub>=0
      print(io," ")
   end
   Printf.@printf(io,"%e  |  ",cs.lb)

   if cs.abs==Inf
      print("          ")
   elseif cs.abs>=0
      print(" ")
   end
   Printf.@printf(io,"%e   ",cs.abs)

   if cs.abs==Inf
      print("          ")
   elseif cs.abs>=0
      print(" ")
   end
   Printf.@printf(io,"%e  |   %e  | %9.3f  %7d",cs.rel,cs.int,cs.time,cs.iter)
end

function InitialConvergenceState()
   ConvergenceState(Inf,-Inf,Inf,Inf,0.0,0,0.0)
end

function has_converged(b::ConvergenceState, a::ConvergenceState)
   if a.time > 1.1*b.time || a.iter > 1.1*b.iter
      true
   elseif a.int <= b.int && (a.abs <= b.abs || a.rel <= b.rel || a.time > b.time || a.iter > b.iter)
      true
   else
      false
   end
end
