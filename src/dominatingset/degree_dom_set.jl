"""
    update_dominated(degree_queue, v, dominated, deleted)

Helper function used to check if a vertex is already dominated and to
update `degree_queue` if it is not.
"""
function update_dominated(
    g::AbstractGraph{T},
    degree_queue::DataStructures.PriorityQueue,
    v::Integer,
    dominated::BitArray{1},
    deleted::BitArray{1}
    ) where T <: Integer

    @inbounds if !dominated[v]
        dominated[v] = true
        if !deleted[v] 
            degree_queue[v] -= 1
        end
        @inbounds @simd for u in neighbors(g, v)
            if !deleted[u] 
                degree_queue[u] -= 1
            end
        end
    end
end

"""
    degree_dominating_set(g)

Greedy Hueristic to solve Minimum Dominating Set.
### Implementation Notes
Initially, all vertices are undominated.
Itertively chooses the vertex that would dominate the most undominated vertices.
### Performance
O( (|V|+|E|)*log(|V|) )
"""
function degree_dominating_set(
    g::AbstractGraph{T}
    ) where T <: Integer 

    nvg = nv(g)  
    dom_set = Vector{T}()  
    deleted = falses(nvg)
    dominated = falses(nvg)
    degree_queue = DataStructures.PriorityQueue(Base.Order.Reverse, zip(collect(1:nv(g)), broadcast(+, degree(g), 1)))

    while !DataStructures.isempty(degree_queue) && DataStructures.peek(degree_queue)[2] > 0
        v = DataStructures.dequeue!(degree_queue)
        deleted[v] = true
        push!(dom_set, v)

        update_dominated(g, degree_queue, v, dominated, deleted)
        for u in neighbors(g, v)
    		update_dominated(g, degree_queue, u, dominated, deleted)
    	end
    end
    return dom_set
end
