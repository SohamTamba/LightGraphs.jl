using DistributedArrays

# Parts of this code were taken / derived from Graphs.jl. See LICENSE for
# licensing details.


"""
    struct FloydWarshallState{T, U}

An [`AbstractPathState`](@ref) designed for Floyd-Warshall shortest-paths calculations.
"""
struct FloydWarshallState{T,U<:Integer} <: AbstractPathState
    dists::Matrix{T}
    parents::Matrix{U}
end

"""
    struct FloydWarshallState{T, U}

A struct used to store data of vertex.
"""
struct vertex_data{T,U<:Integer}
    dists::Vector{T}
    parents::Vector{U}
end

""" 
    LG_FW_it(distr_data, FW_pivot, FW_pivot_data)    

The vertex_data of all nodes will be distributed on workers
Workers will perform an iteration of Floyd Warshall with FW_pivot
as the pivot element for that iteration.
"""
function LG_FW_it(distr_data::DArray{vertex_data{T, U}}, FW_pivot::U, FW_pivot_data::vertex_data{T, U}) where T where U        
    pivot_dists = FW_pivot_data.dists
    pivot_parents = FW_pivot_data.parents
    n_v = size(pivot_dists)[1]
    local_data = distr_data[:L]

    for i in 1:size(local_data)[1]
        dists = local_data[i].dists
        parents = local_data[i].parents
        for v in 1:n_v
            if dists[FW_pivot] == typemax(T) || pivot_dists[v] == typemax(T)
                ans = typemax(T)
            else
                ans = dists[FW_pivot] + pivot_dists[v]
            end
            if dists[v] > ans
                dists[v] = ans
                parents[v] = pivot_parents[v]
            end
        end
        local_data[i] = vertex_data(dists, parents)
    end
    distr_data[:L] = local_data
end


@doc_str """
floyd_warshall_shortest_paths(g, distmx=weights(g))
Use the [Floyd-Warshall algorithm](http://en.wikipedia.org/wiki/Floyd–Warshall_algorithm)
to compute the shortest paths between all pairs of vertices in graph `g` using an
optional distance matrix `distmx`. Return a [`LightGraphs.FloydWarshallState`](@ref) with relevant
traversal information.

### Performance
Space complexity is on the order of ``\\mathcal{O}(|V|^2)``.
"""
function floyd_warshall_shortest_paths(
    g::AbstractGraph{U},
    distmx::AbstractMatrix{T} = weights(g)
) where T where U
    n_v = nv(g)
    dists = fill(typemax(T), (Int(n_v), Int(n_v)))
    parents = zeros(U, (Int(n_v), Int(n_v)))

    for v in 1:n_v
        dists[v, v] = zero(T)
    end
    undirected = !is_directed(g)
    for e in edges(g)
        u = src(e)
        v = dst(e)

        d = distmx[u, v]

        dists[u, v] = min(d, dists[u, v])
        parents[u, v] = u
        if undirected
            dists[v, u] = min(d, dists[v, u])
            parents[v, u] = v
        end
    end

    @eval @everywhere vertex_data = $vertex_data
    @eval @everywhere LG_FW_it = $LG_FW_it
    
    FW_data = [ vertex_data([dists[i, j] for i in 1:n_v], [parents[i, j] for i in 1:n_v]) for j in 1:n_v]
    distr_data = distribute(FW_data, procs=procs())

    #When many workers present, a load balanced broadcast of the remotecall should be used
    for pivot in vertices(g)
        wait_list = Vector{Future}(nworkers())
        pivot_data = distr_data[pivot]
        for p in 1:nworkers()
           wait_list[p] = remotecall(LG_FW_it, p+1, distr_data, pivot, pivot_data)
        end
        LG_FW_it(distr_data, pivot, pivot_data)
        for p in 1:nworkers()
            wait(wait_list[p])
        end
    end

    dists = Matrix{T}(n_v, n_v)
    parents = Matrix{U}(n_v, n_v)

    FW_data = convert(Array{vertex_data}, distr_data) #Gather distributed data into main worker

    for i in 1:n_v
        dists[i, :] = FW_data[i].dists 
        parents[i, :] = FW_data[i].parents
    end

    return FloydWarshallState(dists, parents)
end




function enumerate_paths(s::FloydWarshallState{T,U}, v::Integer) where T where U<:Integer
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


enumerate_paths(s::FloydWarshallState) = [enumerate_paths(s, v) for v in 1:size(s.parents, 1)]
enumerate_paths(st::FloydWarshallState, s::Integer, d::Integer) = enumerate_paths(st, s)[d]
