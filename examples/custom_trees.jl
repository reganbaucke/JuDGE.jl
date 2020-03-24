using JuDGE

leafnodesA = [ [1,1,1,1],
               [1,1,2,1],
               [1,2,1,1,1],
               [1,2,1,1,2],
               [1,2,1,2],
               [1,2,1,3] ]

probsA = [0.2,0.3,0.1,0.2,0.05,0.15]

treeA, probabilityA = tree_from_leaves(leafnodesA,probsA)

nodesA = collect(treeA)

probabilityA(nodesA[10])
probabilityA(get_node(treeA,leafnodesA[3]))
probabilityA(get_node(treeA,[1,2,1,1,1]))
JuDGE.print_tree(treeA)


nodeprobsB = [1.0,
               [0.5,
                 [0.4,
                   [1.0]],
                 [0.6,
                   [1.0]]],
               [0.5,
                 [1.0,
                   [0.6,
                     [1/3],
                     [2/3]],
                   [0.1],
                   [0.3]]]]

treeB,probabilityB = tree_from_nodes(nodeprobsB)

nodesB = collect(treeB)

probabilityB(treeB.children[2].children[1].children[1].children[1])
probabilityB(nodesB[10])
JuDGE.print_tree(treeB)

treeC,probabilityC = tree_from_file(joinpath(@__DIR__,"treeC.tree"))

nodesC = collect(treeC)

probabilityC(get_node(treeC,[1,2,1,1,1]))
JuDGE.print_tree(treeC)
