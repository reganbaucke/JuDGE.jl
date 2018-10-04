@everywhere(include("count_heads.jl"))

a = @spawn count_heads(100)
b = @spawn count_heads(100)
c = 3

function thisfunction(this::Int64)
    this + c
end
