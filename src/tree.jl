
# trait for iterablility
abstract type IterableTrait end
struct Iterable{T} <: IterableTrait end
struct NotIterable{T} <: IterableTrait end
IterableTrait(::Type) = NotIterable{Any}()
IterableTrait(::Type{A}) where {A<:AbstractArray{T}} where {T} = Iterable{T}()
IterableTrait(::Type{A}) where {A<:Dict{U,T}} where {T} where {U} = Iterable{T}()

# definition of a tree
abstract type AbstractTree end

mutable struct Leaf <: AbstractTree
    name::String
    function Leaf()
        return new("")
    end
    function Leaf(name::String)
        return new(name)
    end
end

mutable struct Tree <: AbstractTree
    children::Any
    name::String
    Tree(children) = Tree(IterableTrait(typeof(children)), children)
    function Tree(::Iterable{<:AbstractTree}, children)
        new(children, "")
    end
end

function Base.copy(leaf::Leaf)
    Leaf()
end

function Base.copy(tree::Tree)
    Tree(map(copy, tree.children))
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

function Base.show(io::IO, tree::AbstractTree)
    if typeof(tree) == Tree
        if tree.name != ""
            print(
                io,
                "Subtree rooted at node " *
                tree.name *
                " containing " *
                string(count(tree)) *
                " nodes",
            )
        else
            print(io, "Subtree containing " * string(count(tree)) * " nodes")
        end
    else
        if tree.name != ""
            print(io, "Leaf node " * tree.name)
        else
            print(io, "Leaf node")
        end
    end
end

function count(tree::AbstractTree)
    i = 1
    if typeof(tree) == Tree
        for child in tree.children
            i += count(child)
        end
    end
    return i
end


function print_tree(some_tree)
    function helper(tree::Tree, depth)
        println("  "^depth * "--" * tree.name)
        for child in tree.children
            helper(child, depth + 1)
        end
    end
    function helper(leaf::Leaf, depth)
        println("  "^depth * "--" * leaf.name)
    end
    helper(some_tree, 0)
    nothing
end

function print_tree(some_tree,pr::Dict{AbstractTree,Float64})
    function helper(tree::Tree, depth)
        println("  "^depth * "--" * tree.name * " (" * string(pr[tree]) * ")")
        for child in tree.children
            helper(child, depth + 1)
        end
    end
    function helper(leaf::Leaf, depth)
        println("  "^depth * "--" * leaf.name * " (" * string(pr[leaf]) * ")")
    end
    helper(some_tree, 0)
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
            helper(child, collection)
        end
    end
    result = Array{AbstractTree,1}()
    helper(tree, result)
    result
end

function Base.collect(leaf::Leaf)
    [leaf]
end

function label_nodes(tree::AbstractTree)
    hist = history2(tree)
    function helper(leaf::Leaf, collection)
        leaf.name = hist(leaf)
        push!(collection, leaf)
    end
    function helper(someTree::Tree, collection)
        someTree.name = hist(someTree)
        push!(collection, someTree)
        for child in someTree.children
            helper(child, collection)
        end
    end
    result = Array{AbstractTree,1}()
    helper(tree, result)
    result
end

### collect all the nodes of the tree into an array in a depth first fashion
function get_leafnodes(tree::AbstractTree)
    function helper(leaf::Leaf, collection)
        push!(collection, leaf)
    end
    function helper(someTree::Tree, collection)
        for child in someTree.children
            helper(child, collection)
        end
    end
    result = Array{Leaf,1}()
    helper(tree, result)
    result
end
##
# Given a tree1, this function returns a function which takes in a subtree and returns its parent in the context of the given tree
function parent_builder(tree::T where {T<:AbstractTree})
    function helper(super::Tree, some_tree::AbstractTree)
        if super == some_tree
            return nothing
        end
        for child in super.children
            if child == some_tree
                return super
            end
        end
        just_checking = map(x -> helper(x, some_tree), super.children)
        for result in just_checking
            if result != nothing
                return result
            end
        end
    end
    function helper(super::Leaf, some_tree)
        nothing
    end
    return x -> helper(tree, x)
end


##
# Given a tree, this function returns a function which gives Conditionally Uniform Probabilities values for subtrees of this tree
function ConditionallyUniformProbabilities(tree::T where {T<:AbstractTree})
    parentfunction = parent_builder(tree)
    function result(subtree::T where {T<:AbstractTree})
        if subtree == tree
            return 1.0
        else
            parent = parentfunction(subtree)
            return result(parent) / length(parent.children)
        end
    end
    prob=Dict{AbstractTree,Float64}()
    for node in collect(tree)
        prob[node]=result(node)
    end
    return prob
end

function convert_probabilities(tree::T where {T<:AbstractTree}, probabilities::Dict{AbstractTree,Float64})
    parentfunction = parent_builder(tree)
    function result(subtree::T where {T<:AbstractTree})
        if subtree == tree
            return probabilities[subtree]
        else
            parent = parentfunction(subtree)
            return result(parent) * probabilities[subtree]
        end
    end
    prob=Dict{AbstractTree,Float64}()
    for node in collect(tree)
        prob[node]=result(node)
    end
    return prob
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
function history(tree::T where {T<:AbstractTree})
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
    x -> helper(Array{AbstractTree,1}(), x)
end

function history2(tree::T where {T<:AbstractTree})
    parents = parent_builder(tree)
    function helper(state, subtree)
        if subtree == tree
            state = "1" * state
            return state
        else
            parent = parents(subtree)
            for i = 1:length(parent.children)
                if subtree == parent.children[i]
                    state = string(i) * state
                    break
                end
            end
            state = helper(state, parents(subtree))
        end
        state
    end
    x -> helper("", x)
end

# This function builds a tree from a depth and degree.
function narytree(depth::Int64, degree::Int64)
    function helper(height)
        if height == 0
            output = Leaf()
        else
            v = Array{AbstractTree,1}()
            for i = 1:degree
                push!(v, helper(height - 1))
            end
            output = Tree(v)
        end
        output
    end
    tree=helper(depth)
    label_nodes(tree)
    tree
end

function get_node(tree::AbstractTree, indices::Array{Int64,1})
    node = tree
    for i = 2:length(indices)
        node = node.children[indices[i]]
    end
    return node
end

mutable struct Node
    children::Array{Node,1}
    pr::Float64
    name::String
end

function save_tree_to_file(tree::AbstractTree,filename::String)
    nodes=collect(tree)
    parent=parent_builder(tree)

    f = open(filename,"w")
    println(f,"n,p")
    for n in nodes
        if parent(n)==nothing
            println(f,n.name*","*"-")
        else
            println(f,n.name*","*parent(n).name)
        end
    end
    close(f)
end

# Construct tree from Array of leaf nodes, and corresponding probabilities
function tree_from_leaves(leafnodes::Array{Array{Int64,1},1}, probs::Array{Float64,1})
    prob = Dict{AbstractTree,Float64}()

    function groupnode(node::Node)
        if length(node.children) == 0
            output = Leaf()
        else
            v = Array{AbstractTree,1}()
            for i = 1:length(node.children)
                push!(v, groupnode(node.children[i]))
                node.pr += node.children[i].pr
            end
            output = Tree(v)

        end
        prob[output] = node.pr
        output
    end

    root = Node(Array{Node,1}(), 0.0, "")
    for i = 1:length(leafnodes)
        if leafnodes[i][1] != 1
            error("there must be a unique root node")
        end
        parent = root
        for j = 2:length(leafnodes[i])
            while leafnodes[i][j] > length(parent.children)
                n = Node(Array{Node,1}(), 0.0, "")
                push!(parent.children, n)
            end
            parent = parent.children[leafnodes[i][j]]
            if j == length(leafnodes[i])
                parent.pr = probs[i]
            end
        end
    end
    tree = groupnode(root)
    label_nodes(tree)
    return (tree, prob)
end

# Construct tree from Array of leaf nodes, without probabilities
function tree_from_leaves(leafnodes::Array{Array{Int64,1},1})
    function groupnode(node::Node)
        if length(node.children) == 0
            output = Leaf()
        else
            v = Array{AbstractTree,1}()
            for i = 1:length(node.children)
                push!(v, groupnode(node.children[i]))
            end
            output = Tree(v)

        end
        output
    end

    root = Node(Array{Node,1}(), 0.0, "")
    for i = 1:length(leafnodes)
        if leafnodes[i][1] != 1
            error("there must be a unique root node")
        end
        parent = root
        for j = 2:length(leafnodes[i])
            while leafnodes[i][j] > length(parent.children)
                n = Node(Array{Node,1}(), 0.0, "")
                push!(parent.children, n)
            end
            parent = parent.children[leafnodes[i][j]]
            if j == length(leafnodes[i])
            end
        end
    end

    tree = groupnode(root)
    label_nodes(tree)
    return tree
end

# Construct tree from nested Vector, the first element of vector is the
# probability of a node, remaining elements are vectors representing child nodes
function tree_from_nodes(nodes::Vector{Any})
    prob = Dict{AbstractTree,Float64}()

    function groupnode(node::Any, prob, prev)
        if length(node) == 1
            output = Leaf()
        else
            v = Array{AbstractTree,1}()
            for i = 2:length(node)
                push!(v, groupnode(node[i], prob, prev * node[1]))
            end
            output = Tree(v)
        end
        prob[output] = node[1] * prev
        output
    end
    tree = groupnode(nodes, prob, 1.0)
    label_nodes(tree)
    return (tree, prob)
end

# Construct tree from a file, each line in the file is of the form B,A,...
# representing an arc in the tree, from node "A" to node "B". The total
# number of columns is arbitrary, these columns are converted into dictionaries
# with the key being the node.
function tree_from_file(filename::String)
    data=Dict{Symbol,Dict{Node,Float64}}()
    data2=Dict{Symbol,Dict{AbstractTree,Float64}}()
    nodes=Dict{String,Node}()
    count=Dict{Node,Int64}()
    first=true
    f=open(filename)
    headers=Array{Symbol,1}()
    for l in eachline(f)
        if first
            first=false
            a=split(l,",")
            if a[1]!="n" || a[2]!="p"
                error("First row of tree file should be \"n,p,...\"")
            end
            for i in 3:length(a)
                push!(headers,Symbol(a[i]))
                data[Symbol(a[i])]=Dict{Node,Float64}()
                data2[Symbol(a[i])]=Dict{AbstractTree,Float64}()
            end
        else
            a=split(l,",")
            if length(a)==0
                @warn("Found blank line")
            elseif length(a)!=length(headers)+2
                @warn("Found line with wrong number of elements: \""*l*"\"")
                continue
            end

            if !((string)(a[2]) in keys(nodes))
                if (string)(a[2])!="-"
                    nodes[(string)(a[2])]=Node(Array{Node,1}(),1.0,(string)(a[2]))
                    count[nodes[(string)(a[2])]]=0
                end
            end
            if !((string)(a[1]) in keys(nodes))
                nodes[(string)(a[1])]=Node(Array{Node,1}(),1.0,(string)(a[1]))
                count[nodes[(string)(a[1])]]=0
            end
            if (string)(a[2])!="-"
                push!(nodes[(string)(a[2])].children,nodes[(string)(a[1])])
            end
            for i in 3:length(a)
                data[headers[i-2]][nodes[(string)(a[1])]]=parse(Float64, a[i])
            end
        end
    end

    close(f)

    for (nn,n) in nodes
        for i in n.children
            count[i]+=1
            if count[i]>1
                error("Node "*i.name*" defined multiple times")
            end
        end
    end

    found=0
    root=nothing
    for n in keys(count)
        if count[n]==0
            root=n
            found+=1
        end
    end

    if found==0
        error("No root node found")
    elseif found>1
        error("Multiple root nodes found")
    end

    function groupnode(node::JuDGE.Node)
        if length(node.children)==0
            output=Leaf(node.name)
        else
            v=Array{AbstractTree,1}()
            for n in node.children
                push!(v,groupnode(n))
            end
            output=Tree(v)
            output.name=node.name
        end
        for h in headers
            data2[h][output]=data[h][node]
        end
        output
    end

    tree=groupnode(root)
    return (tree,data2)
end
