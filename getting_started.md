# Getting started


## Requirements

JuDGE requires Julia-0.6.4, and a valid Gurobi installation. For academics,
Gurobi provides a free academic license.


## Installation

JuDGE is installed by the `Pkg` utility provided by Julia. In the Julia REPL,
simply make the following function call.
    Pkg.clone("https://github.com/reganbaucke/JuDGE.jl")

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
of the root node of the tree, and a list of all the nodes in the tree. A node
consists of

- it's parent
- it's children (if any)
- the probability of arriving at the node
- problem data relating to the node

For now, we will build a tree of depth 3, where each node has 2 children with
uniform probabilities using `buildtree`:

    mytree = buildtree(depth=3,degree=2)
`mytree` is now a tree which contains 7 nodes.

Let us now populate our tree with the problem data. Let's say the initial
capacity of our knapsack is 10 units, while the volume of the possible expansion
of the knapsack is 5 units.

    initial_capacity = 10
    volume_of_expansion = 5

We will declare a structure that will help us organize the remaining data for our problem.

    mutable struct Knapsack
        itemreward::Array{Float64,1}
        volume::Array{Float64,1}
        investcost
    end

For our instance of the problem, we will use the following data: for each node,
we will have a knapsack problem with three items to choose from, each with
different rewards, and different volumes.

    # rows is nodes, columns is items
    itemreward = [
        1.0 3.0 8.0; # Node [1] (the root node)
        2.0 4.0 8.0; # Node[1,1]
        3.0 7.0 9.0; # Node[1,2]
        3.0 7.0 9.0; # Node[1,1,1]
        4.0 7.0 7.0; # Node[1,1,2]
        4.0 7.0 9.0; # Node[1,2,1]
        3.0 2.0 1.0 ]# Node[1,2,2]

    # rows is nodes, columns is items
    itemvolume = [
        2.0 4.0 4.0;
        1.0 9.0 1.0;
        2.0 6.0 7.0;
        1.5 7.0 9.0;
        5.5 9.0 9.0;
        6.5 6.0 9.0;
        5.5 6.0 9.0 ]

Finally, the cost of expanding the knapsack at a given point in the tree is
given by:

    # rows is nodes
    investcost = [
        10;
        4;
        10;
        3;
        3;
        3;
        3 ]

We now populate the `data` field for each node in the following manner:

    for n in 1:length(mytree.nodes)
        mytree.nodes[n].data = Knapsack(itemreward[n,:], itemvolume[n,:], investcost[n])
    end

Now with our tree built, and populated with our problem data, we can initialize
an empty new `JuDGEModel` based on our tree:

    m = JuDGEModel(mytree)


We need to provide our `JuDGEModel` with three pieces of information in order to
build and solve correctly. These are the definition of: `JuDGEexpansion!`,
`JuDGEsubproblem!`, and `JuDGEexpansionscosts!`.

Firstly, we need to declare to JuDGE that we are going to have one *expansion
variable* in our problem.

    JuDGEexpansions!(m) do sp
        @expansion(sp,bag_extension)
    end
If you are familiar with the JuMP syntax, you will notice that the `@expansion`
macro mimics that of the `@variable` macro from JuMP. For the time being, JuDGE
only supports binary expansion variables.

Next, we will tell JuDGE the price of expanding the bag at a certain node in the
tree.

    JuDGEexpansioncosts!(m) do master,n,expansion
        # bring expansions into scope, this makes it easier to write the following expression
        bag_extension = expansion[:bag_extension]

        return @expression(master,n.data.investcost*bag_extension)
    end
Here we are telling JuDGE that: in the master problem, price the `bag_extension`
at node `n` at exactly the value stored in `n.data.investcost`. Notice that this
function must return a *JuMP expression*.

Finally, we need to tell JuDGE what optimization problem is occuring at each
node in the `JuDGETree`, we do this by the `JuDGEsubproblem!` function

    JuDGEsubproblem!(m) do sp,n,expansion
        # bring expansions into scope for subproblems, this is useful for writing the JuMP constraints etc.
        bag_extension = expansion[:bag_extension]

        items = 1:5

        # set up the sub problem variables
        @variable(sp, y[items], category=:Bin)

        #set up objective: note that expansions aren't costed here
        @objective(sp, Min, sum(-n.data.itemreward[i]*y[i] for i in items))

        # set up the constraints
        @constraint(sp, sum(n.data.volume[i]*y[i] for i in items) <= initial_capacity + bag_extension*volume_of_expansion )
    end
The optimization problem at each node problem is a classical knapsack problem.
However, notice in the constraint that we are allowing the capacity of our
knapsack to expand by the product of `bag_extension` and `volume_of_extension`.

At this point, we have now constructed a valid `JuDGEModel`.
We can now solve our model by making a call to `JuDGEsolve!`:

    JuDGEsolve!(m,GurobiSolver(OutputFlag=0)) do time, iterations, lb,ub
        if iterations > 200
            return true
        end
        return false
    end
Currently, JuDGE only supports Gurobi as the subproblem and master problem
solvers.

Upon the call of `JuDGEsolve!`, JuDGE will first build the `JuDGEModel`.
Depending on the size of the tree, this can be a time consuming process. Once
built, JuDGE will then solve the `JuDGEModel`.

We can view the optimal solution to our problem by calling

    print_expansions(m)

For more complicated post-solution analysis, we have the ability to inspect
every aspect of our `JuDGEModel` after the call to `JuDGEsolve!`.
