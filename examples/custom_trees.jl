using JuDGE

leafnodes = [ [1,1,1,1],
              [1,1,2,1],
              [1,2,1,1,1],
              [1,2,1,1,2],
              [1,2,1,2],
              [1,2,1,3] ]

probs = [0.2,0.3,0.1,0.2,0.05,0.15]

treeA,probabilityA =tree_from_leaves(leafnodes,probs)

nodesA=collect(treeA)

probabilityA(nodesA[10])
probabilityA(get_node(treeA,leafnodes[3]))
probabilityA(get_node(treeA,[1,2,1,1,1]))
print_tree(treeA)

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

treeB,probabilityB =tree_from_nodes(nodeprobsB)

nodesB=collect(treeB)

probabilityB(treeB.children[2].children[1].children[1].children[1])
probabilityB(nodesB[10])
print_tree(treeB)
