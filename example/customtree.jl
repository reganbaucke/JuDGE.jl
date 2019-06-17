using JuDGE

leafnodes = [ [1,1,1,1],
              [1,1,2,1],
              [1,2,1,1,1],
              [1,2,1,1,2],
              [1,2,1,2],
              [1,2,1,3] ]

probabilities = [0.3,0.1,0.05,0.05,0.1,0.4]
tree = customtree(leafnodes,probabilities)

print(tree)
