
# trait for iterablility
abstract type IterableTrait end
struct Iterable{T} <: IterableTrait end
struct NotIterable{T} <: IterableTrait end
IterableTrait(::Type) = NotIterable{Any}()
IterableTrait(::Type{A}) where {A<:AbstractArray{T}} where {T} = Iterable{T}()
IterableTrait(::Type{A}) where {A<:Dict{U,T}} where {T} where {U} = Iterable{T}()
IterableTrait(::Type{A}) where {A<:T} where {T} = Iterable{T}()


# definition of a tree
abstract type AbstractTree end
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

function Base.map(f, leaf::Leaf)
   f(leaf)
end
function Base.map(f, tree::Tree)
   Tree(map(x -> map(f, x), tree.children))
end

function Base.map(f, dict::Dict)
   Dict(key => f(dict[key]) for key in keys(dict))
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

### collect all the nodes of the tree into an array in a depth first fashion
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
function Base.collect(leaf::Leaf)
   [leaf]
end

##
# Given a tree1, this function returns a function which takes in a subtree and returns its parent in the context of the given tree
function parent_builder(tree::T where {T <: AbstractTree})
   function helper(super::Tree,some_tree::AbstractTree)
      if super == some_tree
         return nothing
      end
      for child in super.children
         if child == some_tree
            return super
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


##
# Given a tree, this function returns a function which gives Conditionally Uniform Probabilities values for subtrees of this tree
function ConditionallyUniformProbabilities(tree::T where {T <: AbstractTree})
   parentfunction = parent_builder(tree)
   function result(subtree::T where T <: AbstractTree)
      if subtree == tree
         return 1.0
      else
         parent = parentfunction(subtree)
         return result(parent)/length(parent.children)
      end
   end
   return result
end

# Given a tree, this function returns a function which takes a subtree and returns its depth in the context of tree
function depth(tree::AbstractTree)
   parents = parent_builder(tree)
   function result(subtree::AbstractTree)
      if tree == subtree
         return 0
      else
         return result(parents(subtree)) + 1
      end
   end
end

# Given a tree, this function returns a function which takes a subtree and returns that subtree's history back up to the root node in the constext of tree
function history(tree::T where T <: AbstractTree)
   parents = parent_builder(tree)
   function helper(state, subtree)
      if subtree == tree
         push!(state, subtree)
         return state
      else
         push!(state, subtree)
         helper(state, parents(subtree))
      end
      state
   end
   x -> helper(Array{AbstractTree,1}(),x)
end

# This function builds a tree from a children generator and a depth.
function narytree(n::Int64,generator)
   f(::Leaf) = Tree(generator())
   out = Leaf()
   for i = 1:n
      out = map(f, out)
   end
   out
end
