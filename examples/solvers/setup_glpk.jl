import GLPK
using JuMP

JuDGE_SP_Solver =
    optimizer_with_attributes(GLPK.Optimizer, "msg_lev" => 0, "mip_gap" => 0.0)

JuDGE_MP_Solver = optimizer_with_attributes(
    # TODO(odow): what on earth is this?
    (method = GLPK.INTERIOR) -> GLPK.Optimizer(),
    "msg_lev" => 0,
    "mip_gap" => 0.0,
)

JuDGE_DE_Solver =
    optimizer_with_attributes(GLPK.Optimizer, "msg_lev" => 2, "mip_gap" => 0.0)
