using SuperFastGraphs
using Base.Test
using LightGraphs

# g = loadgraph("./testdata/barabasi_example.lg", "graph")
g = loadgraph("./testdata/karate.lg", "graph")
g = barabasi_albert(10000, 30)
tic()
t = triangles(g)
toc()

tic()
triangleNumberD = triangleCountingDegree(g)
toc()

println(triangleNumberD == triangles(g))