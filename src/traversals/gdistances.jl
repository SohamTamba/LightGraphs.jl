
"""
    gdistances!(g, source, dists; sort_alg=QuickSort)

Fill `dists` with the geodesic distances of vertices in `g` from source vertex (or
collection of vertices) `source`. `dists` should be a vector of length `nv(g)` 
filled with `typemax(T)`. Return `dists`.

For vertices in disconnected components the default distance is `typemax(T)`.

An optional sorting algorithm may be specified (see Performance section).

### Performance
`gdistances` uses `QuickSort` internally for its default sorting algorithm, since it performs
the best of the algorithms built into Julia Base. However, passing a `RadixSort` (available via
[SortingAlgorithms.jl](https://github.com/JuliaCollections/SortingAlgorithms.jl)) will provide
significant performance improvements on larger graphs.
"""
function gdistances!(g::AbstractGraph{T}, source, vert_level; sort_alg = QuickSort) where T
    n = nv(g)
    visited = zeros(Bool, n)
    n_level = one(T)
    cur_level = Vector{T}()
    sizehint!(cur_level, n)
    next_level = Vector{T}()
    sizehint!(next_level, n)
    @inbounds for s in source
        vert_level[s] = zero(T)
        visited[s] = true
        push!(cur_level, s)
    end
#count_q = count_E = count_b = 0
    while !isempty(cur_level)
        @inbounds for v in cur_level
#count_q+=1
            @inbounds @simd for i in outneighbors(g, v)
#count_E+=1
                if !visited[i]
#count_b+=1
                    push!(next_level, i)
                    vert_level[i] = n_level
                    visited[i] = true
                end
            end
        end
        n_level += one(T)
        empty!(cur_level)
        cur_level, next_level = next_level, cur_level
        sort!(cur_level, alg = sort_alg)
    end
#println("count_q = count_q")
#println("count_E = count_E")
#println("count_b = count_b")
    return vert_level
end

"""
    gdistances(g, source; sort_alg=QuickSort)

Return a vector filled with the geodesic distances of vertices in  `g` from
`source`. If `source` is a collection of vertices each element should be unique.
For vertices in disconnected components the default distance is `typemax(T)`.

An optional sorting algorithm may be specified (see Performance section).

### Performance
`gdistances` uses `QuickSort` internally for its default sorting algorithm, since it performs
the best of the algorithms built into Julia Base. However, passing a `RadixSort` (available via
[SortingAlgorithms.jl](https://github.com/JuliaCollections/SortingAlgorithms.jl)) will provide
significant performance improvements on larger graphs.
"""
gdistances(g::AbstractGraph{T}, source; sort_alg = Base.Sort.QuickSort) where T = gdistances!(g, source, fill(typemax(T), nv(g)); sort_alg = sort_alg)

