abstract type IterableTrait end
struct Iterable{T} <: IterableTrait end
struct NotIterable{T} <: IterableTrait end
IterableTrait(::Type) = NotIterable{Any}()
IterableTrait(::Type{A}) where {A<:AbstractArray{T}} where {T} = Iterable{T}()
IterableTrait(::Type{A}) where {A<:Dict{U,T}} where {T} where {U} = Iterable{T}()
IterableTrait(::Type{A}) where {A<:T} where {T} = Iterable{T}()

const Option{T} = Union{Some{T},Nothing}

# struct Some{T} <: Option{T}
#    value::T
# end


abstract type AbstractTree end

struct IdFunct{T}
   value::T
end
IterableTrait(::Type{A}) where {A<:IdFunct{T}} where {T} = Iterable{T}()

mutable struct Leaf <: AbstractTree end

mutable struct Tree <: AbstractTree
   children
   Tree(children) = Tree(IterableTrait(typeof(children)), children)
   function Tree(::Iterable{<:AbstractTree}, children)
      new(children)
   end
end

function Base.copy(leaf::Leaf)
   Leaf()
end

function Base.copy(tree::Tree)
   Tree(map(copy,tree.children))
end

# tree is a functor
function Base.map(f, leaf::Leaf)
   f(leaf)
end

function Base.map(f, tree::Tree)
   Tree(map(x -> map(f, x), tree.children))
end

function Base.map(f, dict::Dict)
   Dict(key => f(dict[key]) for key in keys(dict))
end

function random_tree()
   function recursive()
      max = random_int(2,6)
      if (random_int(0,10) <= 2)
         Tree(map(x -> recursive(),collect(1:max)))
      else
         Leaf()
      end
   end
   return recursive()
end

function print_tree(some_tree)
   function helper(tree::Tree,depth)
      println("  "^depth * "--")
      for child in tree.children
         helper(child, depth + 1)
      end
   end
   function helper(leaf::Leaf,depth)
      println("  "^depth * "--")
   end
   helper(some_tree,0)
   nothing
end

function random_int(a::Int,b::Int)
   (abs(rand(Int)) % (b-a)) + a
end


function parent_builder(tree::Tree)
   function helper(super::Tree,some_tree::AbstractTree)
      if super == some_tree
         return nothing
      end
      for child in super.children
         if child == some_tree
            return Some(super)
         end
      end
      just_checking = map(x -> helper(x,some_tree), super.children)
      for result in just_checking
         if result != nothing
            return result
         end
      end
   end
   function helper(super::Leaf,some_tree)
      nothing
   end
   return x -> helper(tree,x)
end

c = Leaf()
b = Tree([Leaf(), c, Leaf(), Leaf()])
a = Tree([Leaf(), b])


function Base.collect(tree::Tree)
   function helper(leaf::Leaf, collection)
      push!(collection, leaf)
   end
   function helper(someTree::Tree, collection)
      push!(collection, someTree)
      for child in someTree.children
         helper(child,collection)
      end
   end
   result = Array{AbstractTree,1}()
   helper(tree,result)
   result
end


function ConditionallyUniformProbabilities(tree::AbstractTree)
   parentfunction = parent_builder(tree)
   function result(subtree)
      if subtree == tree
         return 1.0
      else
         parent = parentfunction(subtree).value
         return result(parent)/length(parent.children)
      end
   end
   return result
end

function BespokeProbabilities(tree::T) where T <: AbstractTree
   parentfunction = parent_builder(tree)
   function result(subtree::U) where U <: AbstractTree
      if subtree == tree
         return 1.0
      else
         parent = parentfunction(subtree).value
         if subtree == parent.children[1]
            return result(parent)*0.6666666
         elseif subtree == parent.children[2]
            return result(parent)*0.3333333
         elseif subtree == paren.children[3]
            return result(parent)*0.3333333
         else
            return 0
         end
      end
   end
   return result
end

# write a depth function
function depth(tree::AbstractTree)
   parents = parent_builder(tree)
   function result(subtree::AbstractTree)
      if tree == subtree
         return 0
      else
         return result(parents(subtree).value) + 1
      end
   end
end
somefunction(::Leaf) = Tree([Leaf(), Leaf(), Leaf()])

sometree = somefunction(Leaf())
biggertree = map(somefunction,sometree)
probab = BespokeProbabilities(biggertree)

function builder()
   initial = 10
   function nodal(node)
      model = Model()
      @variable(model, y[1:2])
      @variable(model, u[1:2])
      #things in the right hand side must be the same for every node
      @constraint(model, y[1] + y[2] <= initial + u[1] + u[2])
      @objective(model, Min, y[1] + y[2])
      return model
   end
end

struct JuDGEModel
   tree::AbstractTree
   master_problem::JuMP.Model
   sub_problems::Dict{AbstractTree,JuMP.Model}
   function JuDGEModel(tree,
      probability_function,
      sub_problem_builder)
      this = new()
      tree
      master_problem = JuMP.Model()
      println("got here")
      sub_problems = Dict(i => sub_problem_builder(i) for i in collect(tree))
      return new(tree,master_problem,sub_problems)
   end
end
