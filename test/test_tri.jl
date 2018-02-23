using SuperFastGraphs
using Base.Test
using LightGraphs

# g = loadgraph("./testdata/barabasi_example.lg", "graph")
g = loadgraph("./testdata/karate.lg", "graph")
g = loadgraph("./testdata/internet.lg", "graph")
g = barabasi_albert(1000, 10)

tic()
t = triangles(g)
toc()

tic()
triangleNumberD = triangleCountingDegree(g)
toc()

println(triangleNumberD == triangles(g))