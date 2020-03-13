struct ConvergenceState
   ub::Float64
   lb::Float64
   abs::Float64
   rel::Float64
   time::Float64
   iter::Int64
end

function Base.show(io::IO,cs::ConvergenceState)
      Printf.@printf(io,"%15.3f %15.3f %15.3f %15.3f %10.3f %7d", cs.ub, cs.lb, cs.abs, cs.rel, cs.time, cs.iter)
end

function InitialConvergenceState()
   ConvergenceState(Inf,-Inf,Inf,Inf,0.0,0)
end

function has_converged(b::ConvergenceState, a::ConvergenceState)
   if a.abs <= b.abs
      true
   elseif a.rel <= b.rel
      true
   elseif a.time > b.time
      true
   elseif a.iter > b.iter
      true
   else
      false
   end
end
