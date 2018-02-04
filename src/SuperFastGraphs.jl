module SuperFastGraphs

using LightGraphs
using DataStructures

export sampleDistance
export diameter!
export ccSample
export triangleCounting
export triangleCountingDegree
export prunedBFS

function sampleDistance(g::AbstractGraph)
	numV = nv(g)
	distApp = 0
	k = convert(Int64, floor(log(numV) * 1000)) # 100 per precisione allo 0.1
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
# u : vertex 
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

function prunedBFS(g::AbstractGraph, v::Int64, x::Int64, preProcessData )

	q = Queue(Int64)
	visited = zeros(Int64, nv(g)) #0: not visited, 1:visited, 2:completely explored
	enqueue!(q, v)
	f = Int64[]
	push!(f, v)
	depth = 0
	timeToDepthIncrease = 1
	elementsToDepthIncrease = 0
	nodeAtDepth = Int64[]
	fdValue = 0
	nodeExplored = 0
	ccV = 1
	while v !in preProcessData[ccV] #trovo la componente connessa di v
		ccV++
	end
	while length(q)>0
		e = dequeue!(q)
		push!(nodeAtDepth, e)
		visited[e] = 2 
		neighborhood = neighbors(g, e)
		elementsToDepthIncrease = elementsToDepthIncrease + length(neighborhood)
		timeToDepthIncrease = timeToDepthIncrease - 1
		if timeToDepthIncrease == 0
			#computing bound
			n = len(vertices(g))
			fdValue = fdValue + (length(nodeAtDepth)*depth)
			nodeExplored = nodeExplored + length(nodeAtDepth)
			gamma_d1 = 0
			for node in nodeAtDepth
			    gamma_d1 = gamma_d1 + degree(node)
			end
			bound = max( ((preprocess[2][ccV]-1)^2)/(fdValue - gamma_d1 + (depth +2)*( preprocess[2][ccV]- nodeExplored)) ,
						 ((preprocess[3][ccV]-1)^2)/(fdValue - gamma_d1 + (depth +2)*( preprocess[3][ccV]- nodeExplored))
						) 
						/ (n-1)
			if x > bound
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
	return ((length(preProcessData[ccV]) -1 )^2)/((nv(g)-1)*fdValue)
end

function preProcess(g::AbstractGraph) #preprocess for compute upper boud of closeness centrality
	if (is_directed(g) && is_strongly_connected(g))
		return ([collect(vertices(g))], [nv(g)], [nv(g)])
	elseif (!is_directed(g) && is_connected(g))
		return ([collect(vertices(g))], [nv(g)], [nv(g)])
	elseif (!is_directed(g) && !is_connected(g))
		ccg = connectedComponents(g)
		lenccg = Int64[]
		for cc in ccg
			push!(lenccg, length(cc))
		end
		return (ccg, lenccg, lenccg)
	elseif (is_directed(g) && !is_strongly_connected(g) && is_weakly_connected(g))
		ccg = connectedComponents(g)
		ccDAGTopologicalSorted = topological_sort_by_dfs(condensation(g))
		alpha = zeros(Int64, length(ccDAGTopologicalSorted))
		omega = zeros(Int64, length(ccDAGTopologicalSorted))
		omega[end] = length(ccg[end])
		alpha[end] = length(ccg[end])
		maxAlpha = alpha[end]
		for i in length(ccDAGTopologicalSorted)-1:1 #dynamic programming
			omega[i] = length(ccg[i])+omega[i+1]
			alpha[i] = maxAlpha + length(ccg[i])
			if (alpha[i] > maxAlpha)
				maxAlpha = alpha[i]
			end
		end
		return (ccg, alpha, omega)
	end
end

function fastClosenessCentrality(g::AbstractGraph)
    preprocessData = preProcess(g)
    xk = 0
    topk = Int64[]
    for v in vertices(g)
    	t_xk = prunedBFS(g, v, xk, preProcessData)
    	if t_xk != 0
    		xk = t_xk
    		push!(topk, v)
    	end
    end
    #TODO: order by degree
end

end # module
