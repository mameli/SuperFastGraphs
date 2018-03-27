using LightGraphs
using SuperFastGraphs

g = loadgraph("test/testdata/kite.lg", "kite")
println(topKcc(g, 10))

