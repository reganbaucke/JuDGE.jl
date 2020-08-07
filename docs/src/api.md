# API Reference

## AbstractTree Functions

### Defining Trees
```@docs
JuDGE.narytree
JuDGE.tree_from_file
JuDGE.tree_from_leaves
JuDGE.print_tree
```

### Nodes of Trees
```@docs
Base.collect
JuDGE.get_leafnodes
JuDGE.get_node
```

### Tree Probabilities
```@docs
JuDGE.convert_probabilities
JuDGE.ConditionallyUniformProbabilities
JuDGE.UniformLeafProbabilities
```

### Other Tree functions
```@docs
JuDGE.depth
JuDGE.history
JuDGE.parent_builder
```

## JuDGE Functions
```@docs
JuDGE.JuDGEModel
JuDGE.solve(::JuDGEModel)
JuDGE.branch_and_price
JuDGE.constraint_branch
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

## JuDGE Output
```@docs
JuDGE.write_solution_to_file(::JuDGEModel,::String)
JuDGE.print_expansions(::JuDGEModel)
```

## Deterministic Equivalent
```@docs
JuDGE.DetEqModel
JuDGE.solve(::DetEqModel)
```

## Deterministic Equivalent Output
```@docs
JuDGE.write_solution_to_file(::DetEqModel,::String)
JuDGE.print_expansions(::DetEqModel)
```
