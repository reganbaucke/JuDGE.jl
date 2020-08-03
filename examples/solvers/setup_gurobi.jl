using Gurobi

env = Gurobi.Env()

JuDGE_MP_Solver = optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0,
																		 "Method" => 2,
																		 "Crossover" => 0,
																		 "MIPGap" => 0.0)

JuDGE_SP_Solver = optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0,
																		 "MIPGap" => 0.0)

JuDGE_DE_Solver = optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 1,
 																		 "MIPGap" => 0.0)
