
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
    ID::Union{Nothing,AbstractTree}
    name::String
    parent::Union{Nothing,AbstractTree}
    ext::Dict{Symbol,Any}
    function Leaf()
        return new(nothing,"",nothing,Dict{Symbol,Any}())
    end
    function Leaf(parent::AbstractTree)
        return new(nothing,"",parent,Dict{Symbol,Any}())
    end
    function Leaf(name::String)
        return new(nothing,name,nothing,Dict{Symbol,Any}())
    end
    function Leaf(name::String,parent::AbstractTree)
        return new(nothing,name,parent,Dict{Symbol,Any}())
    end
end

mutable struct Tree <: AbstractTree
    children::Union{Array{AbstractTree,1},Array{Leaf,1},Array{Tree,1}}
    ID::Union{Nothing,AbstractTree}
    name::String
    parent::Union{Nothing,AbstractTree}
    ext::Dict{Symbol,Any}
    Tree(children) = Tree(IterableTrait(typeof(children)), children)
    function Tree(::Iterable{<:AbstractTree}, children)
        new(children, nothing, "",nothing,Dict{Symbol,Any}())
    end
    function Tree()
        new(Array{AbstractTree,1}(),nothing,"",nothing,Dict{Symbol,Any}())
    end
    function Tree(parent::AbstractTree)
        new(Array{AbstractTree,1}(),nothing,"",parent,Dict{Symbol,Any}())
    end
    function Tree(name::String)
        new(Array{AbstractTree,1}(),nothing,name,nothing,Dict{Symbol,Any}())
    end
    function Tree(name::String,parent::AbstractTree)
        new(Array{AbstractTree,1}(),nothing,name,parent,Dict{Symbol,Any}())
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

function getID(node::AbstractTree)
    if node.ID==nothing
        return node
    else
        return node.ID
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
	print_tree(some_tree::AbstractTree, data::Dict{AbstractTree,Any})

Given `some_tree`, this function prints a simple representation of the tree to the REPL.

### Required Arguments
`some_tree` is the tree we wish to visualise

### Optional Arguments
`data` is a dictionary indexed by each node in `some_tree`
"""
function print_tree(some_tree::AbstractTree, data::Dict{AbstractTree,T} where T <: Any)
    function helper(tree::Tree, depth)
        if tree in keys(data)
            println("  "^depth * "--" * tree.name * " (" * string(data[tree]) * ")")
        elseif getID(tree) in keys(data)
            println("  "^depth * "--" * tree.name * " (" * string(data[getID(tree)]) * ")")
        else
            println("  "^depth * "--" * tree.name)
        end

        for child in tree.children
            helper(child, depth + 1)
        end
    end
    function helper(leaf::Leaf, depth)
        if leaf in keys(data)
            println("  "^depth * "--" * leaf.name * " (" * string(data[leaf]) * ")")
        elseif getID(leaf) in keys(data)
            println("  "^depth * "--" * leaf.name * " (" * string(data[getID(leaf)]) * ")")
        else
            println("  "^depth * "--" * leaf.name)
        end
    end
    helper(some_tree, 0)

    nothing
end

"""
	visualize_tree(some_tree::AbstractTree,
        data::Union{Dict{Symbol,Any},Dict{Symbol,Dict{AbstractTree,Float64}}};
        scale_edges=nothing,
        scale_all=1.0)

Given `some_tree`, this function generates a html/js visualization of the tree.

### Required Arguments
`some_tree` is the tree we wish to visualise.

`data` is a dictionary of the data we wish to display, each element is another
dictionary indexed by the nodes of the tree.

### Optional Arguments
`scale_edges` this scales the lengths of the arcs.

`scale_all` this scales the whole network.
"""
function visualize_tree(some_tree::AbstractTree, data::Union{Dict{Symbol,Any},Dict{Symbol,Dict{AbstractTree,Float64}}};scale_edges=nothing,scale_all=1.0)
    maxdata=Dict{Symbol,Any}()
    mindata=Dict{Symbol,Any}()

    function node_json(node)
        temp="{"
        temp*="id:"
        temp*=string(get_id[node])
        temp*=",label:\""
        temp*=node.name*"\""
        temp*=",posX:"
        temp*=string(position[node][1])
        temp*=",posY:"
        temp*=string(position[node][2])
        temp*=",data:{"
        first=true
        for sym in keys(data)
            if typeof(collect(keys(data[sym]))[1]) <: AbstractTree
                if node in keys(data[sym])
                    if first
                        first=false
                    else
                        temp*=","
                    end
                    temp*=string(sym)
                    temp*=":"
                    if typeof(data[sym][node])==String
                        temp*="\""
                        temp*=data[sym][node]
                        temp*="\""
                    else
                        if sym ∉ keys(maxdata)
                            maxdata[sym]=data[sym][node]
                            mindata[sym]=data[sym][node]
                        else
                            if data[sym][node] < mindata[sym]
                                mindata[sym]=data[sym][node]
                            elseif data[sym][node] > maxdata[sym]
                                maxdata[sym]=data[sym][node]
                            end
                        end
                        temp*=string(data[sym][node])
                    end
                end
            elseif typeof(collect(keys(data[sym]))[1]) == String
                first2=true
                for key in keys(data[sym])
                    if node in keys(data[sym][key])
                        if first
                            first=false
                        else
                            temp*=","
                        end
                        if first2
                            first2=false
                            # temp*=","
                            temp*=string(sym)
                            temp*=":{"
                        end
                        temp*=""
                        temp*=replace(key, "," => "_")
                        temp*=":"
                        if typeof(data[sym][key][node])==String
                            temp*="\""
                            temp*=data[sym][key][node]
                            temp*="\""
                        else
                            if sym ∉ keys(maxdata)
                                maxdata[sym]=data[sym][key][node]
                                mindata[sym]=data[sym][key][node]
                            else
                                if data[sym][key][node] < mindata[sym]
                                    mindata[sym]=data[sym][key][node]
                                elseif data[sym][key][node] > maxdata[sym]
                                    maxdata[sym]=data[sym][key][node]
                                end
                            end
                            temp*=string(data[sym][key][node])
                        end
                    end
                end
                if !first2
                    temp*="}"
                end
            end
        end
        temp*="},level:"
        temp*=string(depth(node))
        temp*=",},"
    end

    function arc_json(node,parent)
        temp="{"
        temp*="from:"
        temp*=string(get_id[parent])
        temp*=",to:"
        temp*=string(get_id[node])
        temp*=",},"
    end

    index=1
    nodes="var nodes = ["
    arcs="var edges = ["
    get_id=Dict{AbstractTree,Int64}()

    angles=Dict{AbstractTree,Float64}()
    position=Dict{AbstractTree,Tuple{Float64,Float64}}()

    function setpositions(node::AbstractTree,alpha::Float64,l::Float64,scale::Float64;first=false)
        if typeof(node)==Leaf
            return
        end
        a=2*pi/(length(node.children)-1+2*alpha)

        current=angles[node]+pi+a*alpha
        for child in node.children
            angles[child]=current
            position[child]=(position[node][1]+l*cos(current),position[node][2]+l*sin(current))
            setpositions(child,alpha,l*scale,scale)
            current+=a
        end
    end

    scale_factors = [ [1.0],
                      [1.0,0.9,0.85,0.78,0.73,0.7,0.67,0.65],
                      [1.0,0.65,0.52,0.48],
                      [1.0,0.45,0.42,0.42,0.41,0.4],
                      [1.0,0.44,0.37,0.36],
                      [1.0,0.42,0.34,0.33],
                      [1.0,0.35,0.3],
                      [1.0,0.3,0.26],
                      [1.0,0.27,0.23],
                      [1.0,0.24,0.22] ]

    if scale_edges==nothing
        if typeof(some_tree)==Leaf
            scale_edges=1.0
        else
            dg=length(some_tree.children)
            if dg<=10
                dp=min(depth(collect(some_tree)[end]),length(scale_factors[dg]))
                scale_edges=scale_factors[dg][dp]
            else
                scale_edges=0.22*0.91^(dg-10)
            end
        end
    end

    angles[some_tree]=0.0
    position[some_tree]=(0.0,0.0)

    setpositions(some_tree,0.5,700.0*scale_all,scale_edges,first=true)

    for node in collect(some_tree)
        get_id[node]=index
        nodes*=node_json(node)
        parent=node.parent
        if parent!=nothing
            arcs*=arc_json(node,parent)
        end
        index+=1
    end

    nodes*="];"
    arcs*="];"

    scale="var scale = {"

    for sym in keys(maxdata)
        if length(scale)!=13
            scale*=","
        end
        scale*=string(sym)
        scale*=":"
        scale*="{min:"
        scale*=string(mindata[sym])
        scale*=",max:"
        if maxdata[sym]==mindata[sym]
            maxdata[sym]+=1
        end
        scale*=string(maxdata[sym])
        scale*="}"
    end
    scale*="};"
    min_size=(1.0+79.0/scale_all)*(scale_edges^(0.65/scale_all))^(depth(collect(some_tree)[end]))*scale_all
    max_size=(1.0+79.0/scale_all)*(scale_all)
    data=">"*nodes*"\n"*arcs*"\n"*scale*"\n"*"var node_scale="*string(scale_edges^(0.65/scale_all))*";"*"var min_size="*string(min_size)*";"*"var max_size="*string(max_size)*";"

    s = read(joinpath(dirname(@__DIR__),"visualise","visualize.html"), String)
    c = ">"*read(joinpath(dirname(@__DIR__),"visualise","colors.js"), String)
    filename=joinpath("vis"*string(Int(round((time()*100)%100000)))*".html")
    file=open(filename,"w")
    println(file,replace(replace(s," src=\"smallnetwork.js\">" => data)," src=\"colors.js\">" => c))
    close(file)

    if Sys.iswindows()
        run(`$(ENV["COMSPEC"]) /c start $(filename)`)
    elseif Sys.isapple()
        run(`open $(filename)`)
    elseif Sys.islinux() || Sys.isbsd()
        run(`xdg-open $(filename)`)
    else
        error("Unable to show plot. Try opening the file $(filename) manually.")
    end
    return
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

function Base.collect(leaf::Leaf;order=:depth)
    [leaf]
end

function label_nodes(tree::AbstractTree)
    for t in collect(tree)
        t.name=history2(t)
    end
    tree
end

"""
	get_leafnodes(tree::AbstractTree)

Given `tree`, this function returns an array of corresponding `Leaf` nodes.

### Required Arguments
`tree` is the tree from which we wish to collect leaf nodes

### Example
    leafnodes = JuDGE.get_leafnodes(tree)
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
	ConditionallyUniformProbabilities(tree::AbstractTree)

Given a tree, this function returns a dictionary which maps nodes of the tree to probabilities, given that there are conditionally uniform probabilities over the children of any node.

### Required Arguments
`tree` is the tree for which the probability distribution will be generated

### Example
    probs = ConditionallyUniformProbabilities(tree)
"""
function ConditionallyUniformProbabilities(tree::T where {T<:AbstractTree})
    function result(subtree::T where {T<:AbstractTree})
        if subtree == tree
            return 1.0
        else
            parent = subtree.parent
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
	convert_probabilities(tree::AbstractTree, probabilities::Dict{NodeID,Float64})

Given a dictionary of conditional probabilities for each node in `tree`, this function
returns a dictionary that maps the `NodeID` of each node of `tree` to the corresponding
unconditional probability.

### Required Arguments
`tree` is the tree that the probabilities pertain to

`probabilities` is a dictionary of condition probabilities for each node in `tree`

### Example
    probs = JuDGE.convert_probabilities(tree,probabilities)
"""
function convert_probabilities(tree::T where {T<:AbstractTree}, probabilities::Dict{AbstractTree,Float64})
    function result(subtree::T where {T<:AbstractTree})
        if subtree == tree
            return probabilities[getID(subtree)]
        else
            parent = subtree.parent
            return result(parent) * probabilities[getID(subtree)]
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

Given `tree`, this function returns the depth.
The root node has a depth of 0.

### Required Arguments
`tree` is the node we wish to find the depth for.

### Example
    depth = JuDGE.depth(tree) #returns the depth of a node in a tree
"""
function depth(tree::AbstractTree)
    length(history(tree))-1
end

"""
	history(tree::AbstractTree)

Given `tree`, this function returns the history back up to the root node of `tree`.

### Required Arguments
`tree` is the node that we wish to find the history for.

### Example
    history = JuDGE.history(tree) #get a vector of nodes that precede tree
"""
function history(tree::T where {T<:AbstractTree})
    function helper(state, subtree)
        if subtree.parent == nothing
            push!(state, subtree)
            return state
        else
            push!(state, subtree)
            helper(state, subtree.parent)
        end
        state
    end
    helper(Array{AbstractTree,1}(), tree)
end

function history2(tree::T where {T<:AbstractTree})
    function helper(state, subtree)
        if subtree.parent == nothing
            state = "1" * state
            return state
        else
            parent = subtree.parent
            for i = 1:length(parent.children)
                if subtree == parent.children[i]
                    state = string(i) * state
                    break
                end
            end
            state = helper(state, subtree.parent)
        end
        state
    end
    helper("", tree)
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
    if depth==0
        return Leaf("1")
    end

    tree=Tree("1")
    layer=[tree]
    for dp in 1:depth
        nextlayer=AbstractTree[]
        for t in layer
            for dg in 1:degree
                if dp!=depth
                    nt=Tree(t.name*string(dg),t)
                    push!(t.children,nt)
                    push!(nextlayer,nt)
                else
                    push!(t.children,Leaf(t.name*string(dg),t))
                end
            end
        end
        layer=nextlayer
    end
    #label_nodes(tree)
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

"""
	get_groups(tree::AbstractTree; combine=0)

Given a `tree`, this function will split it up into an array of subtrees. These can
be provided as `blocks` for the `JuDGE.solve()` function to perform partial pricing.

### Required Arguments
`tree` is the tree from which we are finding the node

### Optional Arguments
`combine` this parameter determines the size of the subtrees (higher creates larger subtrees).
If set to 0, the subtrees will be the sets of paths to the leaf nodes.

### Examples
    blocks = get_groups(tree)
"""
function get_groups(tree::AbstractTree; combine=0)
    leafnodes=get_leafnodes(tree)

    groups=Array{AbstractTree,1}[]
    nodes=AbstractTree[]

    for leaf in leafnodes
        node=history(leaf)[combine+1]
        if node ∉ nodes
            push!(nodes,node)
        end
    end

    for node in nodes
        scenario=history(node)
        children = collect(node)
        popfirst!(children)
        scenario = [scenario;children]
        push!(groups,scenario)
    end

    return groups
end

function get_scenarios(tree::AbstractTree)
    scenarios=get_groups(tree,combine=0)

    trees=AbstractTree[]

    for i in 1:length(scenarios)
        t=Leaf()
        t.ID=scenarios[i][1]
        t.name=t.ID.name
        for j in 2:length(scenarios[i])
            t=Tree([t])
            t.ID=scenarios[i][j]
            t.name=t.ID.name
        end
        push!(trees,t)
    end

    return trees
end

mutable struct Node
    children::Array{Node,1}
    pr::Float64
    name::String
end

function save_tree_to_file(tree::AbstractTree,filename::String)
    nodes=collect(tree,order=:breadth)

    f = open(filename,"w")
    println(f,"n,p")
    for n in nodes
        if parent(n)==nothing
            println(f,n.name*","*"-")
        else
            println(f,n.name*","*n.parent.name)
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
            for c in output.children
                c.parent=output
            end
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
            for c in output.children
                c.parent=output
            end
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

# Construct tree from nested Vector, the first element of each vector is the
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
            for c in output.children
                c.parent=output
            end
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

    function groupnode(node::Node)
        if length(node.children)==0
            output=Leaf(node.name)
        else
            v=Array{AbstractTree,1}()
            for n in node.children
                push!(v,groupnode(n))
            end
            output=Tree(v)
            output.name=node.name
            for c in output.children
                c.parent=output
            end
        end
        for h in headers
            data2[h][output]=data[h][node]
        end
        output
    end

    tree=groupnode(root)
    return (tree,data2)
end
