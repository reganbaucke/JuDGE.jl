
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

function print_tree(some_tree::AbstractTree)
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

"""
	print_tree(some_tree::AbstractTree, data::T where T <: Dict)

Given `some_tree`, this function prints a simple representation of the tree to the REPL.

### Required Arguments
`some_tree` is the tree we wish to visualise

### Optional Arguments
`data` is a dictionary indexed by the nodes of `some_tree`
"""
function print_tree(some_tree::AbstractTree, data::T where T <: Dict)
    function helper(tree::Tree, depth)
        println("  "^depth * "--" * tree.name * " (" * string(data[tree]) * ")")
        for child in tree.children
            helper(child, depth + 1)
        end
    end
    function helper(leaf::Leaf, depth)
        println("  "^depth * "--" * leaf.name * " (" * string(data[leaf]) * ")")
    end
    helper(some_tree, 0)
    nothing
end

"""
	collect(tree::Tree;order=:depth)

Given `tree`, this function returns an array of corresponding nodes. By default this will be in a depth-first order.

### Required Arguments
`tree` is the tree from which we wish to collect the nodes

### Optional Arguments
`order` can be set to `:depth` or `:breadth` to specify the order that the nodes are listed in the array.

### Examples
    nodes = collect(tree) #gets an array of nodes from tree in depth-first order
    nodes = collect(tree,order=:breadth) #gets an array of nodes from tree in breadth-first order
"""
function Base.collect(tree::Tree;order=:depth)
    index=1
    collection=Array{AbstractTree,1}()
    push!(collection,tree)
    while index<=length(collection)
        if typeof(collection[index])==Tree
            if order==:depth
                for i in eachindex(collection[index].children)
                    insert!(collection,index+i,collection[index].children[i])
                end
            elseif order==:breadth
                append!(collection,collection[index].children)
            else
                error("\'order\' should be set to :depth or :breadth")
            end
        end
        index+=1
    end
    collection
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

"""
	get_leafnodes(tree::AbstractTree)

Given `tree`, this function returns an array of corresponding `Leaf` nodes.

### Required Arguments
`tree` is the tree from which we wish to collect leaf nodes

### Example
    leafnodes = JuDGE.get_leafnodes(tree) #define a function that returns the parent of each node in the tree
"""
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

"""
	parent_builder(tree::AbstractTree)

Given `tree`, this function returns a function that takes a node of `tree` and returns that node's parent.

### Required Arguments
`tree` is the tree that the parent function will correspond to.

### Example
    parent_fn = JuDGE.parent_builder(tree) #define a function that returns the parent of each node in the tree
    p = parent_fn(node) #get the parent of node
"""
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

"""
	ConditionallyUniformProbabilities(tree::AbstractTree)

Given a tree, this function returns a dictionary which maps nodes of the tree to probabilities, given that there are conditionally uniform probabilities over the children of any node.

### Required Arguments
`tree` is the tree for which the probability distribution will be generated

### Example
    probs = ConditionallyUniformProbabilities(tree)
"""
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

"""
	UniformLeafProbabilities(tree::AbstractTree)

Given a tree, this function returns a dictionary which maps nodes of the tree to probabilities, given that there is a uniform distribution over the leaf nodes

### Required Arguments
`tree` is the tree for which the probability distribution will be generated

### Example
    probs = UniformLeafProbabilities(tree)
"""
function UniformLeafProbabilities(tree::T where {T<:AbstractTree})
    parentfunction = parent_builder(tree)
    function result(subtree::T where {T<:AbstractTree})
        if typeof(subtree) == Leaf
            return p
        else
            pr=0
            for child in subtree.children
                pr+=result(child)
            end
            return pr
        end
    end

    leafnodes=get_leafnodes(tree)
    p=1/length(leafnodes)
    prob=Dict{AbstractTree,Float64}()
    for node in collect(tree)
        prob[node]=result(node)
    end
    return prob
end

"""
	convert_probabilities(tree::AbstractTree, probabilities::Dict{AbstractTree,Float64})

Given a dictionary of conditional probabilities for each node in `tree`, this function returns a dictionary that maps each node of `tree` to the corresponding unconditional probability.

### Required Arguments
`tree` is the tree that the probabilities pertain to

`probabilities` is a dictionary of condition probabilities for each node in `tree`

### Example
    probs = JuDGE.convert_probabilities(tree,probabilities)
"""
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

"""
	depth(tree::AbstractTree)

Given `tree`, this function returns a function that takes a node of `tree` and returns that node's depth.
The root node has a depth of 0.

### Required Arguments
`tree` is the tree that the depth function will correspond to.

### Example
    depth_fn = JuDGE.depth(tree) #define a function that returns the depth of each node in the tree
    dpth = depth_fn(node) #get the depth of node in tree
"""
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

"""
	history(tree::AbstractTree)

Given `tree`, this function returns a function that takes a node of `tree` and returns that node's history back up to the root node of `tree`.

### Required Arguments
`tree` is the tree that the history function will correspond to.

### Example
    history_fn = JuDGE.history(tree) #define a function that returns the history of each node in the tree

    past = history_fn(node) #get a vector of nodes that precede node in tree
"""
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

"""
	narytree(depth::Int64, degree::Int64)

Given the `depth` and `degree`, this function returns an N-ary tree. Note that a depth of 0 return a single `Leaf` node (which is also the root node of the tree).

### Required Arguments
`depth` is the maximum number of arcs from the root node any `Leaf` node

`degree` is the number of children of all nodes, other than the `Leaf` nodes

### Example
    tree = narytree(2,2)
"""
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

"""
	get_node(tree::AbstractTree, indices::Array{Int64,1})

Given a `tree`, and an array of `indices`, this function returns the corresponding node in the tree.

### Required Arguments
`tree` is the tree from which we are finding the node

`indicies` is an array of integer indices identifying a node within `tree`.

### Examples
    node = get_node(tree,[1]) #get the root node
    node = get_node(tree,[1,1]) #get the first child of the root node
    node = get_node(tree,[1,2]) #get the second child of the root node
"""
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

"""
	tree_from_leaves(leafnodes::Array{Array{Int64,1},1}, probs::Array{Float64,1})

Construct tree from Array of leaf nodes, and (optionally) the corresponding probabilities

### Required Arguments
`leafnodes` is an array of arrays defining the set of leaf nodes

### Optional Arguments
`probs` is an array of probabilities for the leaf nodes

### Example
    (tree,prob) = tree_from_leaves([[1,1,1],[1,1,2],[1,2,1],[1,2,2]],[0.25,0.25,0.25,0.25])
    tree = tree_from_leaves([[1,1,1],[1,1,2],[1,2,1],[1,2,2]])
"""#
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

"""
	tree_from_file(filename::String)

Construct tree from a file, each line in the file is of the form B,A,...
representing an arc in the tree, from node "A" to node "B". The total
number of columns is arbitrary. The first row of the file should be
n,p,... these column headers are converted into symbols used to index
the `data`. Each column itself converted into a dictionary, indexed by
the node.

### Required Arguments
`string` is the full path of the file containing the tree

### Example
    tree, data = tree_from_file(joinpath(@__DIR__,"tree.csv"))
"""#
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
