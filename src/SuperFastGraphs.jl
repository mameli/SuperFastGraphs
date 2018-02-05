module SuperFastGraphs

using LightGraphs
using DataStructures

export sampleDistance
export diameter!
export ccSample
export triangleCounting
export triangleCountingDegree
export prunedBFS
export preProcess
export fastClosenessCentrality

function sampleDistance(g::AbstractGraph, precision::Int64 )
	numV = nv(g)
	distApp = 0
	k = convert(Int64, floor(log(numV) * precision)) # precision = 100 per precisione allo 0.1
	k = k < numV ? k : numV # il sottoinsieme preso puo essere al massimo grande quanto quello originale
	vertexSample = randKVector(vertices(g), k)
	
	for v in vertexSample
		distApp = distApp + sum(gdistances(g, v))
	end
	agvApp = distApp / (k * (numV - 1))
	return agvApp
end

# Given a vector the function returns a random subvector with k diffent elements 
function randKVector(v::UnitRange{Int64}, k::Int64)
	kv = zeros(Int64, length(v))
	for i = 1:k
		elem = v[rand(1:end)]
		while kv[elem] != 0
			elem = v[rand(1:end)]
		end
		kv[elem] = elem
	end
	return filter(x -> x != 0, kv)
end

# u : starting vertex
function diameter!(g::AbstractGraph, u::Int64)
		deg = degree_centrality(g, normalize=false)
	# u = findmax(deg)[2]
		dist = gdistances(g, u)
		lb = eccentricity(g, u)
	# println(u)
		ub = 2 * lb
		i = lb
		nodeCounter = 1
		while (ub > lb)
				biu = lb
				for v in vertices(g)
						if (dist[v] == i)
								e = eccentricity(g, v)
				# println(nodeCounter, "  " , v)
				# nodeCounter = nodeCounter + 1
								if (e > biu)
										biu = e
								end
								if (biu == ub)
										break
								end
						end
				end
				if (biu > 2 * (i - 1))
						ub = biu
						lb = biu
				else
						lb = biu
						ub = 2 * (i - 1)
				end
				i = i - 1
		end
		return ub
end

# k : sample size
# u : vertex for cc
function ccSample(g::AbstractGraph, k::Int64, u::Int64)
	ccSample = 0
	tempCC = 0
	vertexSample = randKVector(vertices(g), k)
	n = nv(g)   
	dists = gdistances(g, u)  
	for i in 1:length(vertexSample)
		tempCC = tempCC + (dists[vertexSample[i]] * n) / (k * (n - 1))
	end
	ccSample = 1 / tempCC
	return ccSample
end




function triangleCounting(g::AbstractGraph)
		triangleNumber = 0
		# t = zeros(Int64, nv(g))
		# listT = Set()
		for v in vertices(g)
			vNeighbors = neighbors(g, v)
			for u in vNeighbors
				for w in vNeighbors
					if (u < w && has_edge(g, u, w))
						triangleNumber = triangleNumber + 1
						# t[v] = t[v] + 1
						# t[u] = t[u] + 1
						# t[w] = t[w] + 1
						# # println("triangolo tra ", sort([v,w,u]))
						# push!(listT, sort([v,w,u]))
					end
				end
			end
		end
		return triangleNumber/3
end

function triangleCountingDegree(g::AbstractGraph)
	# listT = triangleCounting(g)
	triangleNumber = 0
	# t = zeros(Int64, nv(g))
	# listT2 = Set()
	for v in vertices(g)
		vNeighbors = neighbors(g, v)
		# println("Vicini di ", v, ":", vNeighbors)
		for u in vNeighbors
			if (degree(g, v) < degree(g, u) || ((degree(g, v) == degree(g, u) && v < u)))		
				for w in vNeighbors
					if (w < u && (degree(g, v) < degree(g, w) || ((degree(g, v) == degree(g, w) && v < w))) && has_edge(g, u, w))
						triangleNumber = triangleNumber + 1
						# t[v] = t[v] + 1
						# t[u] = t[u] + 1
						# t[w] = t[w] + 1
						# println("triangolo tra ", sort([v,w,u]))
						# println("nodo scelto ", v)
						# push!(listT2, sort([v,w,u]))
					end
				end
			end
		end
	end
	# println(t)
	# println(setdiff(listT, listT2))
	# println(triangleNumber)
	return triangleNumber
end

function prunedBFS(g::AbstractGraph, v::Int64, x::Float64, preProcessData )
	q = Queue(Int64)
	visited = zeros(Int64, nv(g)) #0: not visited, 1:visited, 2:completely explored
	enqueue!(q, v)
	f = Int64[]
	push!(f, v)
	depth = 1
	timeToDepthIncrease = 1
	elementsToDepthIncrease = 0
	nodeAtDepth = Int64[]
	fdValue = 0
	nodeExplored = 0
	ccV = 1
	while !in(v, preProcessData[1][ccV]) #trovo la componente connessa di v
		ccV = ccV + 1
	end
	println("---VERTICE ", v,": ")
	while length(q)>0
		e = dequeue!(q)
		push!(nodeAtDepth, e)
		visited[e] = 2 
		neighborhood = neighbors(g, e)
		elementsToDepthIncrease = elementsToDepthIncrease + length(neighborhood)
		timeToDepthIncrease = timeToDepthIncrease - 1
		if timeToDepthIncrease == 0
			n = length(vertices(g))
			fdValue = fdValue + (length(nodeAtDepth)*depth)
			println("fdValue: ",fdValue)
			nodeExplored = nodeExplored + length(nodeAtDepth)
			println("n_d: ", nodeExplored)
			gamma_d1 = 0
			for node in nodeAtDepth
				gamma_d1 = gamma_d1 + degree(g, node)
			end
			println("gamma_d1: ",gamma_d1)
			fdTildeValueAlpha = fdValue - gamma_d1 + (depth+2)*(preProcessData[2][ccV] - nodeExplored)
			println("fdtilde: ",fdTildeValueAlpha)
			fdTildeValueOmega = fdValue - gamma_d1 + (depth+2)*(preProcessData[2][ccV] - nodeExplored)
			a = ((preProcessData[2][ccV]-1)^2) / fdTildeValueAlpha
			b = ((preProcessData[3][ccV]-1)^2) / fdTildeValueOmega
			bound = max( a, b ) / (n-1)
			println("bound: ", bound)
			if x > bound
				println("v: ", v, " bound: ", bound, " NON TOPK")
				return 0 #il nodo non è tra i topk
			end
			#
			timeToDepthIncrease = elementsToDepthIncrease
			elementsToDepthIncrease = 0
			depth = depth + 1
			nodeAtDepth = Int64[]
		end 
		for neighbor in neighborhood
			if visited[neighbor] == 0
				visited[neighbor] = 1
				enqueue!(q, neighbor)
				push!(f,neighbor)
			end
		end
	end

	println("RITORNO VALORE->")
	println("fvalue: ", fdValue)
	println("rv: ", preProcessData[1][ccV])
	d = dijkstra_shortest_paths(g,v).dists
    δ = filter(x->x != typemax(x), d)
	return ((length(preProcessData[1][ccV]) - 1 )^2)/((nv(g)-1)*sum(δ))
end

function preProcess(g::AbstractGraph) #preprocess for compute upper boud of closeness centrality
	if (is_directed(g) && is_strongly_connected(g))
		return ([collect(vertices(g))], [nv(g)], [nv(g)])
	elseif (!is_directed(g) && is_connected(g))
		return ([collect(vertices(g))], [nv(g)], [nv(g)])
	elseif (!is_directed(g) && !is_connected(g))
		ccg = connected_components(g)
		lenccg = Int64[]
		for cc in ccg
			push!(lenccg, length(cc))
		end
		return (ccg, lenccg, lenccg)
	elseif (is_directed(g) && !is_strongly_connected(g) && is_weakly_connected(g))
		ccg = strongly_connected_components(g)
		condensationGraph = condensation(g)
		ccDAGTopologicalSorted = topological_sort_by_dfs(condensationGraph)
		ccgTopologicalSorted = Array{Array{Int64}}(length(ccg))

		alpha = zeros(Int64, length(ccDAGTopologicalSorted))
		omega = zeros(Int64, length(ccDAGTopologicalSorted))
		c = ccg[ccDAGTopologicalSorted[end]]
		omega[end] = length(c)
		alpha[end] = length(c)
		maxAlpha = alpha[end]
		ccgTopologicalSorted[end] = c
		for i in length(ccg)-1:-1:1 #dynamic programming
			c = ccg[ccDAGTopologicalSorted[i]]
			# println(i,": ",c," -> ", length(c))
			ccgTopologicalSorted[i] = c
			omega[i] = length(c)+sum(omega)
			alpha[i] = maxAlpha + length(c)
			if (alpha[i] > maxAlpha)
				maxAlpha = alpha[i]
			end
		end
		return (ccgTopologicalSorted, alpha, omega)
	end
end

function fastClosenessCentrality(g::AbstractGraph)
	preprocessData = preProcess(g)
	println("Dati PreProcess: ", preprocessData)
	xk = 0.0
	topk = []
	for v in vertices(g)
		
		if (degree(g, v) != 0)
			#t_xk = prunedBFS(g, v, xk, preprocessData)
			t_xk = prunedBFS(g, v, 0.0, preprocessData)
			println(v," --- analizzato --->", t_xk)
			if t_xk != 0
				xk = t_xk
				push!(topk, (v, xk))
				println(v," --- trovato top k --->", xk)
			end
		end
	end
	return topk
end

end # module
