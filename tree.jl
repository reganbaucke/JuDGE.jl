workspace();
println("This is a message");

# mutable struct Node
#     data
#     children::Array{Node}
#     parent::Node
#     Node() = new()
# end

mutable struct Tree{T}
    nodes::Array{T,1}
    root::T
    function Tree{T}() where T
        mytree = new()
        mytree.nodes = Array{T,1}()
        root = T()
        mytree.root=root
        root.tree = mytree
        push!(mytree.nodes,root)
        return mytree
    end
end

mutable struct Node
    data
    children::Array{Node,1}
    parent::Node
    tree::Tree
    function Node()
        v = new()
        v.children=Array{Node,1}()
        v.parent=v
        return v
    end
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
        dummy = Node()
        dummy.tree = parent.tree
        push!(parent.tree.nodes,dummy)
        join(parent,dummy)
    end
end

function ischild(x::Node,y::Node)
    return in(y,x.children)
end

function join(parent::Node,child::Node)
    assert(parent.tree==child.tree)
    push!(parent.children,child)
    child.parent =parent
    return nothing
end

function wholetree(tree::Tree,f)
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

tree = Tree{Node}()

narrytree(tree.root,2,3);

wholetree(tree,printdata)
wholetree(tree,mydata)

wholetree(tree,printdata)
wholetree(tree,isleaf)

getparents(tree.root.children[1].children[1])[2]==tree.root.children[1]