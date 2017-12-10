module SuperFastGraphs

using LightGraphs

export sampleDistance

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

end # module
