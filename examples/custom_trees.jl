using JuDGE

function TreeA()
    leafnodesA = [
        [1, 1, 1, 1],
        [1, 1, 2, 1],
        [1, 2, 1, 1, 1],
        [1, 2, 1, 1, 2],
        [1, 2, 1, 2],
        [1, 2, 1, 3],
    ]

    probsA = [0.2, 0.3, 0.1, 0.2, 0.05, 0.15]

    treeA, probabilityA = tree_from_leaves(leafnodesA, probsA)

    JuDGE.print_tree(treeA)

    nodesA = collect(treeA)

    return probabilityA[nodesA[10]] == 0.1 &&
           probabilityA[get_node(treeA, leafnodesA[3])] == 0.1 &&
           probabilityA[get_node(treeA, [1, 2, 1, 1, 1])] == 0.1
end

function TreeB()
    nodeprobsB = [
        1.0,
        [0.5, [0.4, [1.0]], [0.6, [1.0]]],
        [0.5, [1.0, [0.6, [1 / 3], [2 / 3]], [0.1], [0.3]]],
    ]

    treeB, probabilityB = tree_from_nodes(nodeprobsB)

    JuDGE.print_tree(treeB)

    nodesB = collect(treeB)
    return probabilityB[treeB.children[2].children[1].children[1].children[1]] ≈
           0.1 && probabilityB[nodesB[10]] ≈ 0.1
end

function TreeC(; visualise = false)
    treeC, data = tree_from_file(joinpath(@__DIR__, "treeC.tree"))

    nodesC = collect(treeC)
    JuDGE.print_tree(treeC, data[:pr])
    probabilityC = JuDGE.convert_probabilities(treeC, data[:pr])
    JuDGE.print_tree(treeC, probabilityC)

    if visualise
        data[:pr_absolute] = probabilityC
        JuDGE.visualize_tree(treeC, data, scale_edges = 0.5)
    end

    return probabilityC[get_node(treeC, [1, 2, 1, 1, 1])] ≈ 0.1
end

if !isdefined(@__MODULE__, :running_tests) || !running_tests
    TreeA()
    TreeB()
    TreeC(visualise = true)
end
