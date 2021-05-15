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

   if cs.rlx_abs==Inf
      print(io,"          ")
   elseif cs.rlx_abs>=0
      print(io," ")
   end
   Printf.@printf(io,"%e   ",cs.rlx_abs)

   if cs.rlx_rel==Inf || isnan(cs.rlx_rel)
      print(io,"          ")
   elseif cs.rlx_rel>=0
      print(io," ")
   end
   Printf.@printf(io,"%e  |   ",cs.rlx_rel)

   if isnan(cs.int)
      print(io,"         ")
   end

   Printf.@printf(io,"%e  | %9.3f  %7d",cs.int,cs.time,cs.iter)
end

function display(cs::ConvergenceState;relaxation=true)
   print("")
   if cs.obj==Inf
      print("          ")
   elseif cs.obj>=0
      print(" ")
   end
   Printf.@printf(" %e  | ",cs.obj)

   if cs.ub==Inf
      print("          ")
   elseif cs.ub>=0
      print(" ")
   end
   Printf.@printf("%e ",cs.ub)

   if cs.lb==-Inf
      print("         ")
   elseif cs.lb>=0
      print(" ")
   end
   Printf.@printf("%e  |  ",cs.lb)

   if relaxation
      temp1=cs.rlx_abs
      temp2=cs.rlx_rel
   else
      temp1=cs.int_abs
      temp2=cs.int_rel
   end

   if temp1==Inf
      print("          ")
   elseif temp1>=0
      print(" ")
   end
   Printf.@printf("%e   ",temp1)

   if temp2==Inf || isnan(temp2)
      print("          ")
   elseif temp2>=0
      print(" ")
   end
   Printf.@printf("%e  |   ",temp2)

   if isnan(cs.int)
      print("         ")
   end

   Printf.@printf("%e  | %9.3f  %7d",cs.int,cs.time,cs.iter)
   if relaxation
      println("")
   else
      println("*")
   end
end

function InitialConvergenceState()
   ConvergenceState(Inf,Inf,-Inf,Inf,Inf,Inf,Inf,0.0,0,0.0)
end

function has_converged(b::ConvergenceState, a::ConvergenceState)
   if a.int_abs <= b.int_abs || a.int_rel <= b.int_rel || a.rlx_abs <= b.rlx_abs || a.rlx_rel <= b.rlx_rel || a.time > b.time || a.iter >= b.iter
      true
   else
      false
   end
end
