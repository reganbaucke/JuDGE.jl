# Getting started


## Requirements

JuDGE requires Julia-1.3+, and a valid Gurobi installation. For academics,
Gurobi provides a free academic license.


## Installation

JuDGE is installed by the `Pkg` utility provided by Julia. In the Julia REPL,
simply make the following function call.

    ] clone "https://github.com/reganbaucke/JuDGE.jl"
    
Then, in your Julia script, use

    using JuDGE
to import the functions from the JuDGE module into the current namespace.


## First model


### Problem description
For our tutorial, we will consider the following optimization problem: Our goal
is minimize the cost of a stochastic sequence of knapsack problems. We represent
the stochastic process with a discrete scenario tree. At each node in the
scenario tree, we solve a knapsack problem. However at any point in the tree, we
have the ability to expand the capacity of our knapsack at a certain cost. Once
the capacity of our bag has been expanded, we are able use the extra volume for
future knapsack problems from this node of the tree forward.

In this optimization problem, we are trading off against the cost of expanding
our knapsack, versus the ability to fit more into our knapsack. Deciding when to
perform the knapsack expansion is the difficult part of this
optimization problem.


### Solving our problem using JuDGE

First things first, we should bring into our workspace JuDGE, JuMP, and Gurobi using

    using JuDGE
    using JuMP
    using Gurobi
Using these libraries, we will now model and solve our optimization problem.


The lifecycle of a `JuDGEModel` is the following:

1. The definition of a `JuDGETree`;
2. defining the subproblems of the `JuDGEModel`;
3. building the `JuDGEModel`;
3. Solving the `JuDGEModel`.

The user's job is both Steps 1 and 2, while JuDGE will automatically perform
Steps 3 and 4.

A `JuDGETree` can be built in many different ways. A `JuDGETree` simply consists
of the root node of the tree, and a list of all the nodes in the tree. This is
defined as a nested set of subtrees, with the final nodes being leaf nodes. Each
subtree simply defines its children, and there are functions that facilitate the
calculation of its parent and the probability of arriving at the node, and the
data that correspondes to the node, can be referenced through dictionaries.

For now, we will build a tree of depth 3, where each node has 2 children with
uniform probabilities using `buildtree`:

    mytree = narytree(2,2)
`mytree` is now a tree which contains 7 nodes, with depth 2, and degree 2.
(A depth of 0, gives only a single leaf node.)

### Problem data
Let us now associate our tree with the problem data.

For our instance of the problem, we will use the following data: for each node,
we will have a knapsack problem with three items to choose from, each with
different rewards, and different volumes. The structure of this data is arbitrary;
JuDGE just needs to be able to access the relevant data, based on the node being
processed (dictionaries or functions are recommended).

    function invest_cost(node)
      if node == get_node(mytree,[1])
         180.0
      elseif node == get_node(mytree,[1,1])
         50.0
      elseif node == get_node(mytree,[1,2])
         60.0
      elseif node == get_node(mytree,[1,1,1])
         40.0
      elseif node == get_node(mytree,[1,1,2])
         60.0
      elseif node == get_node(mytree,[1,2,1])
         10.0
      elseif node == get_node(mytree,[1,2,2])
         10.0
      end
    end

    function item_volume(node)
      if node == get_node(mytree,[1])
         [6, 2, 1, 1, 1]
      elseif node == get_node(mytree,[1,1])
         [8, 2, 2, 2, 1]
      elseif node == get_node(mytree,[1,2])
         [8, 1, 1, 1, 3]
      elseif node == get_node(mytree,[1,1,1])
         [4, 4, 3, 1, 2]
      elseif node == get_node(mytree,[1,1,2])
         [1, 3, 1, 1, 2]
      elseif node == get_node(mytree,[1,2,1])
         [7, 3, 1, 1, 1]
      elseif node == get_node(mytree,[1,2,2])
         [2, 5, 2, 1, 2]
      end
    end

    function item_reward(node)
      if node == get_node(mytree,[1])
         [60, 20, 10, 15, 10]
      elseif node == get_node(mytree,[1,1])
         [8, 10, 20, 20, 10]
      elseif node == get_node(mytree,[1,2])
         [8, 10, 15, 10, 30]
      elseif node == get_node(mytree,[1,1,1])
         [40, 40, 35, 10, 20]
      elseif node == get_node(mytree,[1,1,2])
         [15, 35, 15, 15, 20]
      elseif node == get_node(mytree,[1,2,1])
         [70, 30, 15, 15, 10]
      elseif node == get_node(mytree,[1,2,2])
         [25, 50, 25, 15, 20]
      end
    end

### Subproblems
We need to now define the subproblems. These are JuMP models with some JuDGE-
specific features. For our knapsack problem:

    function sub_problems(node)
      model = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0))
      @expansion(model, bag)
      @expansioncosts(model, bag*invest_cost(node))
      @variable(model, y[1:5], Bin)
      @expansionconstraint(model, BagExtension, sum(y[i]*item_volume(node)[i] for i in 1:5) <= 3 + 4 * bag)
      @sp_objective(model, sum(-item_reward(node)[i] * y[i] for i in 1:5))
      return model
    end

The three elements of this that make it a JuDGE subproblem are:

`@expansion(model, bag)` This defines the expansion variables, and supports
standard JuMP vectorized variable declaration. These will be binary.
      
`@expansioncosts` This declares an expression for the costs of investment; this
must be linear (an AffExpr).

`@expansionconstraint(...)` Only constraints of this type can reference the
expansion variable. The expansion must be on the right-hand side, and the
constraint must be less than or equal to.
     
The overall optimization problem at each node problem is a classical knapsack
problem. We have hard-coded that the initial volume of the knapsack is 3, and
the investment in the bag increases it by 4.

### Solving JuDGE Model
Now with our tree built and the problem data referenced, we can initialize the
`JuDGEModel` based on our tree, subproblems, and solver.

    judy = JuDGEModel(mytree, ConditionallyUniformProbabilities, sub_problems,
           optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0))

`ConditionallyUniformProbabilities` simply applies a uniform conditional probability
distribution for child nodes. Any function mapping nodes to absolute probabilities
can be used here.

At this point, we have now constructed a valid `JuDGEModel`.
We can now solve our model by making a call to `JuDGE.solve`:

    JuDGE.solve(judy)

There are a number of optional stopping critieria that can be set here:
    abstol, reltol, rlx_abstol, rlx_reltol, duration, iter.

Currently, we recommend using JuDGE with Gurobi as the subproblem and master problem
solvers. Any solvers can be specified, but the master problem must return duals, and
a barrier method is recommended to computational effeciency. The subproblems can be
solved with any method, but currently need to be solved to optimality (bound gap of 0).

We can view the optimal solution to our problem by calling

    println("Objective: "*string(objective_value(judy.master_problem)))
    JuDGE.print_expansions(judy,onlynonzero=false)

Finally, if we want to recover the optimal solutions for the nodes, we must fix the
investments and resolve each subproblem.

    JuDGE.fix_expansions(judy)
    println("Re-solved Objective: " * string(JuDGE.resolve_fixed(judy)))
    
    JuDGE.write_solution_to_file(judy,joinpath(@__DIR__,"knapsack_solution.csv"))

## Deterministic Equivalent
A deterministic equivalent can also be automatically constructed using the following
code:

    deteq = DetEqModel(mytree, ConditionallyUniformProbabilities, sub_problems,
            optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0))
    JuDGE.solve(deteq)
    
    println("Deterministic Equivalent Objective: " * string(objective_value(deteq.problem)))
    JuDGE.write_solution_to_file(judy,joinpath(@__DIR__,"knapsack_solution.csv"))
