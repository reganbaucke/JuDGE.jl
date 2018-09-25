# workspace();

# mutable struct Node
#     data
#     children::Array{Node}
#     parent::Node
#     Node() = new()
# end

mutable struct Node{T}
    data
    children::Array{Node,1}
    parent::Node
    tree::T
    function Node{T}(mytree::T = Tree() ) where T
        v = new()
        v.children = Array{Node,1}()
        v.parent = v
        v.tree = mytree
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
    if depth==0
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

function wholetree(f,tree::Tree)
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

function Base.display(tree::Tree)
    println("Tree with $(length(tree.nodes)) nodes")
end

function Base.print(tree::Tree)
    println("Tree with $(length(tree.nodes)) nodes")
end

function Base.show(tree::Tree)
    println("Tree with $(length(tree.nodes)) nodes")
end

function Base.print(node::Node)
    println("Node")
    parents = getparents(node)
    push!(parents,node)
    indices = Int[]
    for n in parents
        push!(indices,getindexofnode(n))
    end
    println(indices)
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

function buildtree(f::Function, depth, splits)
    mytree = Tree()
    # create the root node
    mytree.root = Node(mytree)
    narrytree(mytree.root,depth,splits)
    for n in mytree.nodes
        f(n)
    end
    return mytree
end

function buildtree(depth::Number, splits::Number)
    mytree = Tree()
    # create the root node
    mytree.root = Node(mytree)
    narrytree(mytree.root,depth,splits)
    return mytree
end

# function lineartree
#     array = linear
#     return 
# end

function getnode(tree::Tree,indices::Array{Int64,1})
    node=tree.root
    for i = 2:length(indices)
        node = node.children[indices[i]]
    end
    return node
end

# function getindex(tree::Tree,indices::Array{Int64,1})
#     node=tree.root
#     for i = 2:length(indices)
#         node = node.children[indices[i]]
#     end
#     return node
# end

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

function capacityvariable(var)
    var.meta[:capvar] = :true
end

function capacityconstraint(con)
    con.meta[:capcon] = :true
end

function onlyoneinvestment(con)
    con.meta[:oneinvest] = :true
end