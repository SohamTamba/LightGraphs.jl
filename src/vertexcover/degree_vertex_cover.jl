"""
    degree_vertex_cover(g)

Greedy Hueristic to solve Minimum Vertex Cover.
### Implementation Notes
Itertively chooses the vertex with the largest degree into the covering.
### Performance
O( (|V|+|E|)*log(|V|) )
"""
function degree_vertex_cover(
    g::AbstractGraph{T}
    ) where T <: Integer 

    nvg::T = nv(g)  
    cover = Vector{T}()  
    deleted = falses(nvg)
    degree_queue = DataStructures.PriorityQueue(Base.Order.Reverse, zip(collect(1:nv(g)), degree(g)))

    while !isempty(degree_queue) && DataStructures.peek(degree_queue)[2] > 0
        v = DataStructures.dequeue!(degree_queue)

        deleted[v] = true
        push!(cover, v)
        @inbounds @simd for u in neighbors(g, v)
            if !deleted[u] 
                degree_queue[u] -= 1
            end
        end
    end

    return cover
end
