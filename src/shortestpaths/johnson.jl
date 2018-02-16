"""
    struct JohnsonState{T, U}

An [`AbstractPathState`](@ref) designed for Johnson shortest-paths calculations.
"""

struct JohnsonState{T<:Real, U<:Integer} <: AbstractPathState
    dists::Matrix{T}
    parents::Matrix{U}
end



@doc_str """
johnson_shortest_paths(g, distmx=weights(g); parallel = false)

Use the [Johnson algorithm](https://en.wikipedia.org/wiki/Johnson%27s_algorithm)
to compute the shortest paths between all pairs of vertices in graph `g` using an
optional distance matrix `distmx`.

If the parameter parallel is set true, Dijkstra will run in parallel 
Return a [`LightGraphs.JohnsonState`](@ref) with relevant
traversal information.

Returns the typemin for each dist and 0 for each parent if there is a negative cycle

### Performance
Space complexity is on the order of ``\\mathcal{O}(|V|^2)``.
Assuming Dijkstra is using a min heap for priority queue,
Time complexity is on the order of ``\\mathcal{O}(|V|*|E|*log|V|)``
Parallel implementation will run in ``\\mathcal{O}(|E|*log|V| + |V|*|E|)``

### Dependencies from LightGraphs
bellman_ford_shortest_paths
dijkstra_shortest_paths
parallel_multisource_dijkstra_shortest_paths
"""

function johnson_shortest_paths(
    g::AbstractGraph{U},
    distmx::AbstractMatrix{T} = weights(g);
    parallel::Bool = false
) where T<:Real where U<:Integer

    nvg = nv(g)
    wt_transform = bellman_ford_shortest_paths(g, vertices(g), distmx).dists    

    #Transform weights
    for e in edges(g)
    	distmx[src(e), dst(e)] -= wt_transform[dst(e)] - wt_transform[src(e)]
    end

    dists = Matrix{T}(Int(nvg), Int(nvg))
    parents = Matrix{U}(Int(nvg), Int(nvg))
    
    if !parallel
    	for v in vertices(g)
    		dijk_state = dijkstra_shortest_paths(g, v, distmx)
    		dists[v, :] = dijk_state.dists
      		parents[v, :] = dijk_state.parents
    	end
    else
    	dijk_state = parallel_multisource_dijkstra_shortest_paths(g, vertices(g), distmx)
    	dists = dijk_state.dists
    	parents = dijk_state.parents
    end

    #Undo weight transformation
    for e in edges(g)
    	distmx[src(e), dst(e)] += wt_transform[dst(e)] - wt_transform[src(e)]
    end

    for u in vertices(g), v in vertices(g)
    	dists[u, v] += wt_transform[v] - wt_transform[u]
    end

    result = JohnsonState(dists, parents)
    return result
end

function enumerate_paths(s::JohnsonState{T, U}, v::Integer) where T<:Real where U<:Integer
    pathinfo = s.parents[v, :]
    paths = Vector{Vector{U}}()
    for i in 1:length(pathinfo)
        if (i == v) || (s.dists[v, i] == typemax(T))
            push!(paths, Vector{U}())
        else
            path = Vector{U}()
            currpathindex = i
            while currpathindex != 0
                push!(path, currpathindex)
                currpathindex = pathinfo[currpathindex]
            end
            push!(paths, reverse(path))
        end
    end
    return paths
end

enumerate_paths(s::JohnsonState) = [enumerate_paths(s, v) for v in 1:size(s.parents, 1)]
enumerate_paths(st::JohnsonState, s::Integer, d::Integer) = enumerate_paths(st, s)[d]