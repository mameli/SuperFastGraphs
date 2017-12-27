using SuperFastGraphs
using Base.Test
using LightGraphs

g = loadgraph("./testdata/karate.lg", "graph")

@test sampleDistance(g) == 2.408199643493761

@test diameter!(g, 1) == 5

cc1 = closeness_centrality(g)[1]
cc1Sample = ccSample(g, 10, 1)
@test (cc1Sample < (cc1 + 1)) & (cc1Sample > (cc1 - 1))