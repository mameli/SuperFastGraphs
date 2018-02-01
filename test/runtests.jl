using SuperFastGraphs
using Base.Test
using LightGraphs

# g = Graph(5)

# add_edge!(g, 1, 2)
# add_edge!(g, 1, 3)
# add_edge!(g, 1, 4)
# add_edge!(g, 1, 5)
# # add_edge!(g, 1, 6)
# # add_edge!(g, 1, 7)

# add_edge!(g, 2, 3)
# add_edge!(g, 2, 4)
# # add_edge!(g, 5, 8)
# # add_edge!(g, 2, 9)
# # add_edge!(g, 4, 5)
# # add_edge!(g, 8, 9)

g = loadgraph("./testdata/karate.lg", "graph")

@test sampleDistance(g) == 2.408199643493761

@test diameter!(g, 1) == 5

cc1 = closeness_centrality(g)[1]
cc1Sample = ccSample(g, 10, 1)
@test (cc1Sample < (cc1 + 1)) & (cc1Sample > (cc1 - 1))

triangleNumber = triangleCountingDegree(g)
@test (triangleNumber == (sum(triangles(g))/3))