module SuperFastGraphs

using LightGraphs

export sampleDistance
export diameter!
export ccSample
export triangleCounting
export triangleCountingDegree

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
    	return triangleNumber
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


end # module
