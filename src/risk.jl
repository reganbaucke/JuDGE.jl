struct Risk
	λ::Float64
	α::Float64
	offset::Union{Dict{Leaf,Float64},Nothing}
	bound::Union{Float64,Nothing}
	penalty::Float64
end

function Risk(λ::Float64, α::Float64; offset::Union{Dict{Leaf,Float64},Nothing}=nothing, bound::Union{Float64,Nothing}=nothing, penalty::Float64=100000.0)
	if α<=0.0 || α>1.0
		error("α must be >0 and <=1")
	elseif λ<0.0 || λ>1.0
		error("λ must be >=0 and <=1")
	end
	Risk(λ,α,offset,bound,penalty)
end

function Risk(α::Float64; offset::Union{Dict{Leaf,Float64},Nothing}=nothing, bound::Union{Float64,Nothing}=nothing, penalty::Float64=100000.0)
	if α<=0.0 || α>1.0
		error("α must be >0 and <=1")
	end
	Risk(0.0,α,offset,bound,penalty)
end
