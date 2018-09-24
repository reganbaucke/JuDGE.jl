workspace();

type Node{T}
    data
    tree::T
end

mutable struct Tree
    nodes::Array{Node,1}
    root::Node
    function Tree()
        mytree = new()
        mytree.nodes = Array{Node,1}()
        return mytree
    end
end

function Node(data = nothing, tree::Tree=Tree())
    Node{Tree}(data,tree)
end


