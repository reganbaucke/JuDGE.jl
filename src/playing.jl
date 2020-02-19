abstract type IterableTrait end
struct Iterable{T} <: IterableTrait end
struct NotIterable{T} <: IterableTrait end
IterableTrait(::Type) = NotIterable{Any}()
IterableTrait(::Type{A}) where {A<:AbstractArray{T}} where {T} = Iterable{T}()
IterableTrait(::Type{A}) where {A<:Dict{U,T}} where {T} where {U} = Iterable{T}()
IterableTrait(::Type{A}) where {A<:T} where {T} = Iterable{T}()

abstract type AbstractTree{T} end

struct IdFunct{T}
   value::T
end
IterableTrait(::Type{A}) where {A<:IdFunct{T}} where {T} = Iterable{T}()

mutable struct Leaf{T} <: AbstractTree{T}
   data::T
end

mutable struct Tree{T} <: AbstractTree{T}
   data::T
   children
   Tree(a, children) = Tree(a, IterableTrait(typeof(children)), children)
   function Tree(a::T,::Iterable{<:AbstractTree{T}}, children) where T
      out = new{T}(a,children)
   end
end

function fold(f,init,tree::T) where {T <: AbstractTree}
   function helper(accum, aleaf::Leaf)
      f(accum,aleaf)
   end
   function helper(accum, atree::Tree)
      accum = f(accum, atree)
      for i in atree.children
         accum = helper(accum, i)
      end
      accum
   end
   helper(init,tree)
end

# function fold(f,init,tree::T) where {T <: AbstractTree}
#    function helper(accum, aleaf::Leaf)
#       f(accum,aleaf.data)
#    end
#    function helper(accum, atree::Tree)
#       accum = f(accum, atree.data)
#       for i in atree.children
#          accum = helper(accum, i)
#       end
#       accum
#    end
#    helper(init,tree)
# end

function parent(tree::Tree)
   dic = fold(Dict(),tree) do accum, subtree
      if subtree isa Tree
         for child in subtree.children
            accum[child] = subtree
         end
         accum
      else
         accum
      end
   end
   return (x -> dic[x])
end

function Base.size(tree::T) where T <: AbstractTree
   fold((x,y)-> x + 1, 0, tree)
end

# this function builds a probability function that returns probabilities
# supposing that each split is uniform in probability
function ConditionallyUniform(tree::Tree)
   dic = fold(Dict(),tree) do accum, subtree
      accum[tree] = 1.0
      if subtree isa Tree
         for child in subtree.children
            accum[child] = accum[subtree]/length(subtree.children)
         end
         accum
      else
         accum
      end
   end
   return (x -> dic[x])
end

# this function builds a probability function that returns probabilities
# supposing that each leaf has uniform probability
function LeavesAreUniform(tree::Tree)
end

a = Tree(1,[Leaf(2),Tree(3,[Leaf(4),Leaf(5)])])

# function Base.map(f, leaf::Leaf)
#    Leaf(f(leaf.data))
# end
#
# function Base.map(f, tree::Tree)
#    Tree(f(tree.data),map(x -> map(f, x),tree.children))
# end

function Base.map(f, leaf::Leaf)
   Leaf(f(leaf))
end

function Base.map(f, tree::Tree)
   Tree(f(tree),map(x -> map(f, x),tree.children))
end

function Base.copy(leaf::Leaf)
   Leaf(leaf.data)
end

function Base.copy(tree::Tree)
   Tree(tree.data, map(copy,tree.children))
end

Just(1)


const Option{T} = Union{Some{T},Nothing}

struct Some{T} <: Option{T}
   value::T
end
