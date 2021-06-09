using Test

const _EXAMPLES_DIR = joinpath(dirname(@__DIR__), "examples")

# Setup and initialize the subproblem solvers for the example using GLPK.
running_tests=true
include(joinpath(_EXAMPLES_DIR, "solvers", "setup_glpk.jl"))

@testset "JuDGE tests" begin
	@testset "Multistage newsvendor" begin
		include(joinpath(_EXAMPLES_DIR,"newsvendor.jl"))
		@test newsvendor(cost=5.0,price=8.0,demands=[10,80,100]) ≈ -53.333 atol = 1e-3
		@test newsvendor(cost=5.0,price=8.0,demands=[10,80,100],CVaR=Risk(0.5,0.5)) ≈ -30.0 atol = 1e-3
		@test newsvendor(depth=2,cost=5.0,price=8.0,demands=[10,20,30],CVaR=Risk(0.15,0.5)) ≈ -61.011 atol = 1e-3
		@test newsvendor(depth=3,cost=5.0,price=8.0,demands=[10,20,30],CVaR=Risk(0.15,0.05)) ≈ -90.526 atol = 1e-3
	end

	@testset "Inventory" begin
		include(joinpath(_EXAMPLES_DIR,"inventory.jl"))
		@test inventory(depth=2,degree=2,
						price_array=[0.1172393013694979, 0.25653400961083994, 4.2365322616699785e-6, 0.7161790267880648,
						0.05823720128592225, 0.04993686809222453, 0.9201443039152302]) ≈  -13.949 atol = 1e-3

		@test inventory(depth=2,degree=2,
						price_array=[0.1172393013694979, 0.25653400961083994, 4.2365322616699785e-6, 0.7161790267880648,
						0.05823720128592225, 0.04993686809222453, 0.9201443039152302], risk=Risk(0.1,bound=1.0)) ≈  -10.467 atol = 1e-3
	end

	@testset "Stochastic Knapsack" begin
		include(joinpath(_EXAMPLES_DIR,"knapsacks.jl"))
		@test knapsack_fixed() ≈ -131.25 atol = 1e-3
		@test knapsack_shutdown() ≈ -145.25 atol = 1e-3
		@test knapsack_budget() ≈ -159.0 atol = 1e-3
	end
end

running_tests=false
