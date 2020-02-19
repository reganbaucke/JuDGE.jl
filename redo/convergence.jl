
struct ConvergenceState
   abs::Float64
   rel::Float64
   time::Float64
   iter::Int64
end

function InitialConvergenceState()
   ConvergenceState(Inf,Inf,0.0,0)
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
