using SuperFastGraphs
using Base.Test
using LightGraphs

gd = DiGraph(6)

add_edge!(gd, 1, 2)
add_edge!(gd, 1, 3)
add_edge!(gd, 1, 4)
add_edge!(gd, 5, 1)
add_edge!(gd, 2, 3)
add_edge!(gd, 3, 6)
add_edge!(gd, 4, 5)
#println(preProcess(g))
#
#println("-----------------")
#g = DiGraph(3)
#add_edge!(g, 1, 2)
#add_edge!(g, 2, 3)
#add_edge!(g, 3, 1)
#println(preProcess(g))
#
#println("---------------")
#
# g = Graph(6)
# add_edge!(g, 1, 2)
# add_edge!(g, 1, 3)
# add_edge!(g, 2, 4)
# add_edge!(g, 2, 5)
# add_edge!(g, 5, 4)
# println("connesso: ",is_connected(g))
# println("componenti connesse: ", connected_components(g))
# println(preProcess(g))

println("-------------------")
g = loadgraph("testdata/karate.lg", "graph")
println(g)
println(fastClosenessCentrality(g, 5))
println("-----------------------")
println(closeness_centrality(g))
#println(connected_components(g))