using SuperFastGraphs
using Base.Test
using LightGraphs

# g = loadgraph("./testdata/barabasi_example.lg", "graph")
g = loadgraph("./testdata/karate.lg", "graph")
g = loadgraph("./testdata/internet.lg", "graph")
# g = barabasi_albert(1000, 100)


# g = Graph(5)

# add_edge!(g, 1, 2)
# add_edge!(g, 1, 3)
# add_edge!(g, 1, 4)
# add_edge!(g, 1, 5)
# add_edge!(g, 2, 4)
# add_edge!(g, 2, 3)

tic()
t = triangles(g)
toc()

tic()
triangleNumberD = triangleCountingDegree(g)
toc()

tic()
nt = triangleForward(g)
toc()

println(nt == sum(t)/3)