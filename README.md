# JuDGE

JuDGE stands for: Julia Decomposition for General Expansion. Functionally,
it is a solver which leverages the syntax of the JuMP modelling language to
solve a particular class of optimisation problems.

## Problem Class / Decomposition

JuDGE solves multi-stage stochastic integer programming problems using 
Dantzig-Wolfe decomposition. The user must specify a tree that represents
the uncertainty of the problem, and at each node define a subproblem that
can be a linear or integer program. Further, the expansion variables which
link the subproblems must be declared.

JuDGE automatically generates a master problem and performs column generation
to converge to an optimal solution.

## Stochastic Knapsack Example

JuDGE is distributed with an example of a multi-stage stochastic integer
programming problem. This is a stochastic knapsack problem with investment.
The file knapsack.jl contains the implementation of this problem within the
JuDGE framework, and generateData.jl can be used to create a new (randomized)
data set for the values and sizes of the items in each of the nodes of the
tree.

## Limitations

- For each expansion variable, there can only be one expansion (i.e. variable
is binary).

- JuDGE does not support decisions to reduce capacity.

- The JuDGE master problem presumes natural integrality.

## Bugs

Please raise an `issue' if you experience an error while using JuDGE.
