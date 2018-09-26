module Parallel

using LightGraphs
using LightGraphs: sample, AbstractPathState, JohnsonState, BellmanFordState, FloydWarshallState
using Distributed: @distributed
using Base.Threads: @threads, nthreads, Atomic, atomic_add!, atomic_cas!
using SharedArrays: SharedMatrix, SharedVector, sdata
using ArnoldiMethod
using Random:shuffle
import SparseArrays: sparse
import Base: push!, popfirst!, isempty, getindex

include("shortestpaths/bellman-ford.jl")
include("shortestpaths/dijkstra.jl")
include("shortestpaths/floyd-warshall.jl")
include("shortestpaths/johnson.jl")
include("centrality/betweenness.jl")
include("centrality/closeness.jl")
include("centrality/pagerank.jl")
include("centrality/radiality.jl")
include("centrality/stress.jl")
include("distance.jl")
include("traversals/bfs.jl")
include("traversals/gdistances.jl")
include("traversals/bit_gdistances.jl")
include("traversals/greedy_color.jl")
include("utils.jl")
include("vertexcover/random_vertex_cover.jl")

# Overload until https://github.com/JuliaLang/julia/pull/28651
import Distributed: splitrange
splitrange(N::Integer, np::Integer) = splitrange(Int.(N), np)

end
