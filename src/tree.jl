mutable struct Node{T}
    data
    children::Array{Node,1}
    parent::Node
    tree::T
    p::Float64
    function Node{T}(mytree::T = Tree() ) where T
        v = new()
        v.children = Array{Node,1}()
        v.parent = v
        v.tree = mytree
        v.p = 1;
        push!(mytree.nodes,v)
        return v
    end
end

mutable struct Tree
    nodes::Array{Node,1}
    root::Node
    function Tree()
        mytree = new()
        mytree.nodes = Array{Node,1}()
        return mytree
    end
    function Tree(someNode::Node)
        mytree = new()
        mytree.root = someNode
        mytree.nodes = Array{Node,1}()
        push!(mytree.nodes,someNode)
        return mytree
    end
end

function Node(mytree::Tree = Tree() )
    v = Node{Tree}(mytree)
end

function narrytree(root::Node,depth::Core.Integer,n::Core.Integer)
    if depth==1
        return nothing
    end
    spawn(root,n)
    for child in root.children
        narrytree(child,depth-1,n)
    end
end

function spawn(parent::Node,n::Core.Integer)
    assert(n>0)
    for i = 1:n
        dummy = Node(parent.tree)
        join(parent,dummy)
        dummy.p = parent.p/n
    end
end

function ischild(x::Node,y::Node)
    return in(y,x.children)
end

function join(parent::Node,child::Node)
    assert(parent.tree==child.tree)
    push!(parent.children,child)
    child.parent = parent
    return nothing
end

function wholetree!(f,tree::Tree)
    for node in tree.nodes
        f(node);
    end
end

function mydata(node::Node)
    node.data = 1
    return nothing
end

function printdata(node::Node)
    println(node.data)
    return nothing
end

function isleaf(node::Node)
    println(length(node.children)==0)
end

function getparents(node::Node)
    list = Array{Node,1}();
    while node.parent!=node
        unshift!(list,node.parent)
        node=node.parent
    end
    return list
end

function getindexofnode(node::Node)
    if node.parent == node
        return 1
    end
    for idx in eachindex(node.parent.children)
        if node.parent.children[idx]==node
            return idx
        end
    end
end

function customtree(leafnodes::Array{Array{Int64,1},1},probabilities::Array{Float64,1})
    tree = Tree()
    tree.root=Node(tree)
    for i in 1:length(leafnodes)
        if leafnodes[i][1]!=1
            error("there must be a unique root node")
        end
        parent=tree.root
        for j in 2:length(leafnodes[i])
            while leafnodes[i][j]>length(parent.children)
                n=Node(tree)
                n.parent=parent
                push!(parent.children,n)
            end
            parent=parent.children[leafnodes[i][j]]
        end

        parent.p=probabilities[i]
    end

    determineprobs!(tree.root)

    if tree.root.p != 1
        println("WARNING: root node probability is " * string(tree.root.p))
    end

    return tree
end

function determineprobs!(node)
    if length(node.children)!=0
        temp=0
        for n in node.children
            temp+=determineprobs!(n)
        end
        node.p=temp
    end
    return node.p
end

function Base.display(tree::Tree)
    println("Tree with $(length(tree.nodes)) nodes")
end

function Base.print(tree::Tree)
    println("Tree with $(length(tree.nodes)) nodes")
    for n in tree.nodes
        print(n)
    end
end

function Base.show(tree::Tree)
    println("Tree with $(length(tree.nodes)) nodes")
end

function Base.print(node::Node)
    #println("Node")
    parents = getparents(node)
    push!(parents,node)
    indices = Int[]
    for n in parents
        push!(indices,getindexofnode(n))
    end
    println(string(indices) * ": " * string(node.p))
end

function Base.display(node::Node)
    println("Node")
    parents = getparents(node)
    push!(parents,node)
    indices = Int[]
    for n in parents
        push!(indices,getindexofnode(n))
    end
    println(indices)
end

function Base.show(io::IO,node::Node)
    print(io,"Node")
    parents = getparents(node)
    push!(parents,node)
    indices = Int[]
    for n in parents
        push!(indices,getindexofnode(n))
    end
    print(io,indices)
end

function buildtree(f::Function; depth::Int=0 , degree::Int=0)
    mytree = Tree()
    # create the root node
    mytree.root = Node(mytree)
    narrytree(mytree.root,depth,degree)
    for n in mytree.nodes
        f(n)
    end
    return mytree
end

function buildtree(;depth::Int=0 , degree::Int=0)
    mytree = Tree()
    # create the root node
    mytree.root = Node(mytree)
    narrytree(mytree.root,depth,degree)
    return mytree
end

function Base.start(n::Node{Tree})
    return false
end
Base.next(n::Node{Tree},state)  = (n,true)
Base.done(n::Node{Tree},state)  = state

function getnode(tree::Tree,indices::Array{Int64,1})
    node=tree.root
    for i = 2:length(indices)
        node = node.children[indices[i]]
    end
    return node
end

function lineartree(tree::Tree)
    list = Array{Node,1}()
    push!(list,tree.root)
    for node in  tree.root.children
        push!(list,tree.root)
    end
    for node in  tree.root.children
        push!(list,tree.root)
    end
end

function Base.getindex(tree::Tree,indices...)
    tmp = Array{Int64,1}()
    for i = 1:length(indices)
        push!(tmp, indices[i])
    end
    return getnode(tree,tmp)
end

function stage(node::Node{Tree})
    return length(getparents(node)) + 1
end

function P(n::Node)
    list = Node[]
    list = getparents(n)
    push!(list,n)
end
