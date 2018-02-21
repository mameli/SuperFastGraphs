# SuperFastGraphs

[![Build Status](https://travis-ci.org/mameli/SuperFastGraphs.svg?branch=master)](https://travis-ci.org/mameli/SuperFastGraphs) [![Coverage Status](https://coveralls.io/repos/github/mameli/SuperFastGraphs/badge.svg?branch=master)](https://coveralls.io/github/mameli/SuperFastGraphs?branch=master) [![codecov.io](http://codecov.io/github/mameli/SuperFastGraphs/coverage.svg?branch=master)](http://codecov.io/github/mameli/SuperFastGraphs?branch=master)

## Installing the Package
```
> julia
> Pkg.clone("https://github.com/mameli/SuperFastGraphs.git") 
> Pkg.add("LightGraphs") #  check dependency
> Pkg.add("DataStructures") # check dependency
```

## Usage
```julia
using SuperFastGraphs
using LightGraphs

g = loadgraph("PATH/graph_name.lg", "graph") # create a graph with lightgraphs
```

### Approximate Distance
```julia
sampleDistance(g, 100) 
```

### Faster diameter
```julia
iFub(g) # the starting node will the node with max degree

iFub(g, 4) # the starting node is 4, the execution time can vary with different nodes
```

### Approximate Closeness centrality
```julia
ccSample(g, 10, 1) # calculate the closeness centrality of the vertex 1 with 10 random vertex as sample
```

### Top K Closeness centrality
```julia
topKcc(g, 5) # return the list of 5 vertex with the highest closesess centrality
```
|                     | LightGraphs           | SuperFastGraphs  |
| ------------------  |:---------------------:| ----------------:|
| Distance            | :heavy_check_mark:    |:heavy_check_mark:|
| Distance with prob  | :x:                   |:heavy_check_mark:|
| Diameter            | :heavy_check_mark:    |:heavy_check_mark:|
| ifub                | :x:                   |:heavy_check_mark:|
| Closeness Centrality| :heavy_check_mark:    |:heavy_check_mark:|
| CC with Sampling    | :x:                   |:heavy_check_mark:|
| Top k CC	          | :x:                   |:heavy_check_mark:|
| Triangles           | :heavy_check_mark:    |:heavy_check_mark:|
