# JuDGE API

## JuDGE Functions
```@docs
JuDGE.JuDGEModel
JuDGE.solve(::JuDGEModel)
JuDGE.resolve_subproblems
```

## JuDGE macros for subproblems
```@docs
JuDGE.@expansion
JuDGE.@shutdown
JuDGE.@expansioncosts
JuDGE.@maintenancecosts
JuDGE.@sp_objective
```
## Deterministic Equivalent
```@docs
JuDGE.DetEqModel
JuDGE.solve(::DetEqModel)
```

## Solution Output Functions
### JuDGE Output
```@docs
JuDGE.write_solution_to_file(::JuDGEModel,::String)
JuDGE.print_expansions(::JuDGEModel)
```

### Deterministic Equivalent Output
```@docs
JuDGE.write_solution_to_file(::DetEqModel,::String)
JuDGE.print_expansions(::DetEqModel)
```
