using SuperFastGraphs
using Base.Test
using LightGraphs

gDWeakly = DiGraph(6)

add_edge!(gDWeakly, 1, 2)
add_edge!(gDWeakly, 1, 3)
add_edge!(gDWeakly, 1, 4)
add_edge!(gDWeakly, 5, 1)
add_edge!(gDWeakly, 2, 3)
add_edge!(gDWeakly, 3, 6)
add_edge!(gDWeakly, 4, 5)
println(is_directed(gDWeakly) && !is_strongly_connected(gDWeakly) && is_weakly_connected(gDWeakly))
println(preProcess(gDWeakly))

println("-----------------")
gDStrongly = DiGraph(3)
add_edge!(gDStrongly, 1, 2)
add_edge!(gDStrongly, 2, 3)
add_edge!(gDStrongly, 3, 1)
println(is_directed(gDStrongly) && is_strongly_connected(gDStrongly))
println(preProcess(gDStrongly))

println("---------------")

gNotD_NotC = Graph(6)
add_edge!(gNotD_NotC, 1, 2)
add_edge!(gNotD_NotC, 1, 3)
add_edge!(gNotD_NotC, 2, 4)
add_edge!(gNotD_NotC, 2, 5)
add_edge!(gNotD_NotC, 5, 4)
println("connesso: ",is_connected(gNotD_NotC))
println("componenti connesse: ", connected_components(gNotD_NotC))
println(!is_directed(gNotD_NotC) && !is_connected(gNotD_NotC))
println(preProcess(gNotD_NotC))

println("-------------------")
# g = loadgraph("testdata/karate.lg", "graph")
g = loadgraph("testdata/internet.lg", "graph")
println(g)
tic()
topKcc(g, 33)
toc()
# println(fastClosenessCentrality(g, 5)[1])
println("-----------------------")
tic()
closeness_centrality(g)
toc()
# println(sort!(closeness_centrality(g))[end-5])
#println(connected_components(g))