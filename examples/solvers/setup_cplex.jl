import CPLEX
using JuMP

JuDGE_MP_Solver = optimizer_with_attributes(
    CPLEX.Optimizer,
    "CPXPARAM_SolutionType" => 2,
    "CPXPARAM_MIP_Tolerances_AbsMIPGap" => 0.0,
    "CPXPARAM_MIP_Tolerances_Integrality" => 10^-6,
    "CPX_PARAM_SCRIND" => 0,
    "CPXPARAM_MIP_Display" => 0,
    "CPXPARAM_ParamDisplay" => 0,
)

JuDGE_SP_Solver = optimizer_with_attributes(
    CPLEX.Optimizer,
    "CPXPARAM_MIP_Tolerances_MIPGap" => 0.0,
    "CPXPARAM_MIP_Tolerances_Integrality" => 10^-6,
    "CPX_PARAM_SCRIND" => 0,
    "CPXPARAM_ParamDisplay" => 0,
)

JuDGE_DE_Solver = optimizer_with_attributes(
    CPLEX.Optimizer,
    "CPXPARAM_MIP_Tolerances_MIPGap" => 0.0,
    "CPXPARAM_MIP_Tolerances_Integrality" => 10^-6,
    "CPX_PARAM_SCRIND" => 1,
)
