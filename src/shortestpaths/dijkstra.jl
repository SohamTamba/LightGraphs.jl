"""
    struct DijkstraState{T, U}

An [`AbstractPathState`](@ref) designed for Dijkstra shortest-paths calculations.
"""
struct DijkstraState{T <: Real,U <: Integer} <: AbstractPathState
    parents::Vector{U}
    dists::Vector{T}
    predecessors::Vector{Vector{U}}
    pathcounts::Vector{UInt64}
    closest_vertices::Vector{U}
end

"""
    dijkstra_shortest_paths(g, srcs, distmx=weights(g));

Perform [Dijkstra's algorithm](http://en.wikipedia.org/wiki/Dijkstra%27s_algorithm)
on a graph, computing shortest distances between `srcs` and all other vertices.
Return a [`LightGraphs.DijkstraState`](@ref) that contains various traversal information.

### Optional Arguments
- `allpaths=false`: If true, returns a [`LightGraphs.DijkstraState`](@ref) that keeps track of all
predecessors of a given vertex.

### Performance
If using a sparse matrix for `distmx`, you *may* achieve better performance by passing in a transpose of its sparse transpose.
That is, assuming `D` is the sparse distance matrix:
```
D = transpose(sparse(transpose(D)))
```
Be aware that realizing the sparse transpose of `D` incurs a heavy one-time penalty, so this strategy should only be used
when multiple calls to `dijkstra_shortest_paths` with the distance matrix are planned.

# Examples
```jldoctest
julia> using LightGraphs

julia> ds = dijkstra_shortest_paths(CycleGraph(5), 2);

julia> ds.dists
5-element Array{Int64,1}:
 1
 0
 1
 2
 2

julia> ds = dijkstra_shortest_paths(PathGraph(5), 2);

julia> ds.dists
5-element Array{Int64,1}:
 1
 0
 1
 2
 3
```
"""
function dijkstra_shortest_paths(g::AbstractGraph,
    srcs::Vector{U},
    distmx::AbstractMatrix{T}=weights(g);
    allpaths=false,
    trackvertices=false
    ) where T <: Real where U <: Integer

    nvg = nv(g)
    dists = fill(typemax(T), nvg)
    parents = zeros(U, nvg)
    visited = zeros(Bool, nvg)

    pathcounts = zeros(UInt64, nvg)
    preds = fill(Vector{U}(), nvg)
    H = PriorityQueue{U,T}()
    # fill creates only one array.

    for src in srcs
        dists[src] = zero(T)
        visited[src] = true
        pathcounts[src] = 1
        H[src] = zero(T)
    end

    closest_vertices = Vector{U}()  # Maintains vertices in order of distances from source
    sizehint!(closest_vertices, nvg)

    while !isempty(H)
        u = dequeue!(H)

        if trackvertices
            push!(closest_vertices, u)
        end

        d = dists[u] # Cannot be typemax if `u` is in the queue
        for v in outneighbors(g, u)
            alt = d + distmx[u, v]

            if !visited[v]
                visited[v] = true
                dists[v] = alt
                parents[v] = u

                pathcounts[v] += pathcounts[u]
                if allpaths
                    preds[v] = [u;]
                end
                H[v] = alt
            elseif alt < dists[v]
                dists[v] = alt
                parents[v] = u
                #615
                pathcounts[v] = pathcounts[u]
                if allpaths
                    resize!(preds[v], 1)
                    preds[v][1] = u
                end
                H[v] = alt
            elseif alt == dists[v]
                pathcounts[v] += pathcounts[u]
                if allpaths
                    push!(preds[v], u)
                end
            end
        end
    end

    if trackvertices
        for s in vertices(g)
            if !visited[s]
                push!(closest_vertices, s)
            end
        end
    end

    for src in srcs
        pathcounts[src] = 1
        parents[src] = 0
        empty!(preds[src])
    end

    return DijkstraState{T,U}(parents, dists, preds, pathcounts, closest_vertices)
end

dijkstra_shortest_paths(g::AbstractGraph, src::Integer, distmx::AbstractMatrix=weights(g); allpaths=false, trackvertices=false) =
dijkstra_shortest_paths(g, [src;], distmx; allpaths=allpaths, trackvertices=trackvertices)


export int_dijkstra_shortest_paths

using LightGraphs.IPQ

function int_dijkstra_shortest_paths(g::AbstractGraph,
    srcs::Vector{U},
    distmx::AbstractMatrix{T}=weights(g);
    allpaths=false,
    trackvertices=false
    ) where T <: Real where U <: Integer

    nvg = nv(g)
    dists = fill(typemax(T), nvg)
    parents = zeros(U, nvg)
    visited = zeros(Bool, nvg)

    pathcounts = zeros(UInt64, nvg)
    preds = fill(Vector{U}(), nvg)
    H = Int_PriorityQueue{U,T}(nvg)
    # fill creates only one array.

    for src in srcs
        dists[src] = zero(T)
        visited[src] = true
        pathcounts[src] = 1
        H[src] = zero(T)
    end

    closest_vertices = Vector{U}()  # Maintains vertices in order of distances from source
    sizehint!(closest_vertices, nvg)

    while !isempty(H)
        u = dequeue!(H)

        if trackvertices
            push!(closest_vertices, u)
        end

        d = dists[u] # Cannot be typemax if `u` is in the queue
        for v in outneighbors(g, u)
            alt = d + distmx[u, v]

            if !visited[v]
                visited[v] = true
                dists[v] = alt
                parents[v] = u

                pathcounts[v] += pathcounts[u]
                if allpaths
                    preds[v] = [u;]
                end
                H[v] = alt
            elseif alt < dists[v]
                dists[v] = alt
                parents[v] = u
                #615
                pathcounts[v] = pathcounts[u]
                if allpaths
                    resize!(preds[v], 1)
                    preds[v][1] = u
                end
                H[v] = alt
            elseif alt == dists[v]
                pathcounts[v] += pathcounts[u]
                if allpaths
                    push!(preds[v], u)
                end
            end
        end
    end

    if trackvertices
        for s in vertices(g)
            if !visited[s]
                push!(closest_vertices, s)
            end
        end
    end

    for src in srcs
        pathcounts[src] = 1
        parents[src] = 0
        empty!(preds[src])
    end

    return DijkstraState{T,U}(parents, dists, preds, pathcounts, closest_vertices)
end

int_dijkstra_shortest_paths(g::AbstractGraph, src::Integer, distmx::AbstractMatrix=weights(g); allpaths=false, trackvertices=false) =
int_dijkstra_shortest_paths(g, [src;], distmx; allpaths=allpaths, trackvertices=trackvertices)

export heap_dijkstra_shortest_paths

function heap_dijkstra_shortest_paths(g::AbstractGraph,
    srcs::Vector{U},
    distmx::AbstractMatrix{T}=weights(g);
    allpaths=false,
    trackvertices=false
    ) where T <: Real where U <: Integer

    nvg = nv(g)
    @boundscheck checkbounds(distmx, Base.OneTo(nvg), Base.OneTo(nvg))

    dists = fill(typemax(T), nvg)
    parents = zeros(U, nvg)
    visited = zeros(Bool, nvg)

    pathcounts = zeros(UInt64, nvg)
    preds = fill(Vector{U}(), nvg) #fill creates only one array.
    pq = BinaryMinHeap{Tuple{T, U}}() #heaps avoid costly decrease_key operation

    for src in srcs
        dists[src] = zero(T)
        pathcounts[src] = 1
        push!(pq, (zero(T), src))
        preds[src] = []
    end

    closest_vertices = Vector{U}()  #maintains vertices in order of distances from source
    sizehint!(closest_vertices, nvg)

    @inbounds while !isempty(pq)
        d, u = pop!(pq) #(distance of u, u)
        visited[u] && continue #without using decrease_key, it is possible for the pq to contain multiple copies of u
        visited[u] = true

        if trackvertices
            push!(closest_vertices, u)
        end

        for v in outneighbors(g, u)
#            (distmx[u,v] < 0) && error("Negative edges detected. Use bellman_ford_shortest_paths instead.")
            alt = d + distmx[u, v] #new distance to v via u

            if alt < dists[v] #not visited or found a shorter path to visit
                dists[v] = alt
                parents[v] = u
                pathcounts[v] = pathcounts[u]
                push!(pq, (alt, convert(U,v)))

                if allpaths
                    if isassigned(preds[v], 1)
                        resize!(preds[v], 1)
                        preds[v][1] = u
                    else #first time accessing vertex v
                        preds[v] = [u]
                    end
                end
            elseif alt == dists[v]
                pathcounts[v] += pathcounts[u]
                if allpaths
                    push!(preds[v], u)
                end
            end
        end
    end

    if trackvertices
        @inbounds for s in vertices(g)
            if !visited[s]
                push!(closest_vertices, s)
            end
        end
    end

    return DijkstraState{T,U}(parents, dists, preds, pathcounts, closest_vertices)
end

heap_dijkstra_shortest_paths(g::AbstractGraph, src::Integer, distmx::AbstractMatrix=weights(g); allpaths=false, trackvertices=false) =
heap_dijkstra_shortest_paths(g, [src;], distmx; allpaths=allpaths, trackvertices=trackvertices)
