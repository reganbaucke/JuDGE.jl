struct ConvergenceState
   obj::Float64
   ub::Float64
   lb::Float64
   int_abs::Float64
   int_rel::Float64
   rlx_abs::Float64
   rlx_rel::Float64
   time::Float64
   iter::Int64
   int::Float64
   function ConvergenceState(obj, ub, lb, int_abs, int_rel, rlx_abs, rlx_rel, time, iter, frac)
      return new(obj, ub, lb, int_abs, int_rel, rlx_abs, rlx_rel, time, iter, frac)
   end
   function ConvergenceState(obj, ub, lb, time, iter, frac)
      return new(obj, ub, lb, ub-lb, (ub-lb)/abs(lb), obj-lb, (obj-lb)/abs(lb), time, iter, frac)
   end
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

   if cs.int_abs==Inf
      print(io,"          ")
   elseif cs.int_abs>=0
      print(io," ")
   end
   Printf.@printf(io,"%e   ",cs.int_abs)

   if cs.int_rel==Inf || isnan(cs.int_rel)
      print(io,"          ")
   elseif cs.int_rel>=0
      print(io," ")
   end
   Printf.@printf(io,"%e  |   ",cs.int_rel)

   if isnan(cs.int)
      print(io,"         ")
   end

   Printf.@printf("%e  | %9.3f  %7d",cs.int,cs.time,cs.iter)
end

function InitialConvergenceState()
   ConvergenceState(Inf,Inf,-Inf,Inf,Inf,Inf,Inf,0.0,0,0.0)
end

function has_converged(b::ConvergenceState, a::ConvergenceState)
   if a.int_abs <= b.int_abs || a.int_rel <= b.int_rel || a.rlx_abs <= b.rlx_abs || a.rlx_rel <= b.rlx_rel || a.time > b.time || a.iter > b.iter
      true
   else
      false
   end
end
