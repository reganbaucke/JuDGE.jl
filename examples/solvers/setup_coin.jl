import Clp
import Cbc
using JuMP

JuDGE_SP_Solver = optimizer_with_attributes(
    Cbc.Optimizer,
    "logLevel" => 0,
    "ratioGap" => 0.0,
    "presolve" => "off",
    "preprocess" => "off",
)

JuDGE_MP_Solver = (
    optimizer_with_attributes(Clp.Optimizer, "LogLevel" => 0, "Algorithm" => 4),
    JuDGE_SP_Solver,
)

JuDGE_DE_Solver =
    optimizer_with_attributes(Cbc.Optimizer, "logLevel" => 1, "ratioGap" => 0.0)
