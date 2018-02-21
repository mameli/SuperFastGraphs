using SuperFastGraphs
using Base.Test
using LightGraphs

g = loadgraph("./testdata/internet.lg", "graph")

tic()
t = triangles(g)
toc()
tic()
triangleNumber = triangleCounting(g)
toc()
tic()
triangleNumberD = triangleCountingDegree(g)
toc()
println(triangleNumber == (sum(triangles(g))/3))
println(triangleNumberD == (sum(triangles(g))/3))