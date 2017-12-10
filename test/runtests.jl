using SuperFastGraphs
using Base.Test
using LightGraphs

g = loadgraph("./testdata/karate.lg", "graph")

@test sampleDistance(g) == 2.408199643493761
