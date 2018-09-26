"""
    partition_sources!(queue_list, sources)

Partition `sources` using [`LightGraphs.unweighted_contiguous_partition`](@ref) and place
the i^{th} partition  into `queue_list[i]` and set to empty_list[i] to true if the 
i^{th} partition is empty.
"""
function bit_partition_sources!(
    queue_list::Vector{Vector{T}},
    sources::Vector{<:Integer},
    empty_list::Vector{Bool}
    ) where T<:Integer
    
    partitions = LightGraphs.unweighted_contiguous_partition(length(sources), length(queue_list))
    for (i, p) in enumerate(partitions)
        append!(queue_list[i], sources[p])
        empty_list[i] = isempty(p)
    end
end

export bit_gdistances!, bit_gdistances

"""
    bit_gdistances!(g, sources, vert_level; queue_segment_size=20)
    bit_gdistances!(g, source, vert_level; queue_segment_size=20)

Parallel implementation of [`LightGraphs.gdistances!`](@ref) with dynamic load balancing.

### Optional Arguments
- `queue_segment_size = 20`: It is the number of vertices a thread can claim from a queue
at a time. For graphs with uniform degree, a larger value of `queue_segment_size` could
improve performance.

### References
- [Avoiding Locks and Atomic Instructions in Shared-Memory Parallel BFS Using Optimistic 
Parallelization](https://www.computer.org/csdl/proceedings/ipdpsw/2013/4979/00/4979b628-abs.html).
"""
function bit_gdistances!(
    g::AbstractGraph{T}, 
    sources::Vector{<:Integer},
    vert_level::Vector{T};
    queue_segment_size::Integer=20
    ) where T <:Integer
 
    nvg = nv(g)
    n_t = nthreads()
    segment_size = convert(T, queue_segment_size) # Type stability
    fill!(vert_level, typemax(T))
    visited = zeros(Bool, nvg)
    visited_bit = falses(nvg)

    #bitVector not thread safe
    next_level_t = [sizehint!(Vector{T}(), cld(nvg, n_t)) for _ in Base.OneTo(n_t)]
    cur_level_t = [sizehint!(Vector{T}(), cld(nvg, n_t)) for _ in Base.OneTo(n_t)]
    cur_front_t = ones(T, n_t)
    queue_explored_t = zeros(Bool, n_t)

    for s in sources    
        visited[s] = true
        visited_bit[s] = true
        vert_level[s] = zero(T)
    end
    bit_partition_sources!(cur_level_t, sources, queue_explored_t)
    is_cur_level_t_empty = isempty(sources)
    n_level = zero(T)

#=
count_q = zeros(Int64, n_t)
count_V = zeros(Int64, n_t)
count_E = zeros(Int64, n_t)
count_b = zeros(Int64, n_t)
count_bit = zeros(Int64, n_t)
=#
    while !is_cur_level_t_empty
        n_level += one(T)

        let n_level=n_level # let block used due to bug #15276
        @threads for thread_id in Base.OneTo(n_t)
            #Explore current level in parallel
            @inbounds next_level = next_level_t[thread_id]

            @inbounds for t_range in (thread_id:n_t, 1:(thread_id-1)), t in t_range
                queue_explored_t[t] && continue
                cur_level = cur_level_t[t]
                cur_len = length(cur_level)

                # Explore cur_level_t[t] one segment at a time.
                while true 
                    local_front = cur_front_t[t]  # Data race, but first read always succeeds
                    cur_front_t[t] += segment_size # Failure of increment is acceptable

                    (local_front > cur_len || local_front <= zero(T)) && break                    
                    while local_front <= cur_len && cur_level[local_front] != zero(T)
                        v = cur_level[local_front]
                        cur_level[local_front] = zero(T)
                        local_front += one(T)
#count_q[thread_id]+=one(Int64)
                        # Check if v was successfully read.
                        (visited[v] && vert_level[v] == n_level-one(T)) || continue
#count_V[thread_id]+=one(Int64)
                        for i in outneighbors(g, v)
#count_E[thread_id]+=one(Int64)
                        if !visited_bit[i]
                            visited_bit[i] = true
#count_bit[thread_id]+=one(Int64)
                            # Data race, but first read on visited[i] always succeeds
                            if !visited[i]
                                visited[i] = true
#count_b[thread_id]+=one(Int64)
                                vert_level[i] = n_level
                                #Concurrent visited[i] = true always succeeds
                                push!(next_level, i)
                            end
                        end
                        end
                    end   
                end
                queue_explored_t[t] = true        
            end      
        end
        end

        is_cur_level_t_empty = true
        @inbounds for t in Base.OneTo(n_t)
            cur_level_t[t], next_level_t[t] = next_level_t[t], cur_level_t[t]
            cur_front_t[t] = one(T)
            empty!(next_level_t[t])
            queue_explored_t[t] = isempty(cur_level_t[t])
            is_cur_level_t_empty = is_cur_level_t_empty && queue_explored_t[t]
        end               
    end
#=
println("count_q = $(count_q)")
println(sum(count_q))
println("count_V = $(count_V)")
println(sum(count_V))
println("count_E = $(count_E)")
println(sum(count_E))
println("count_bit = $(count_bit)")
println(sum(count_bit))
println("count_b = $(count_b)")
println(sum(count_b))
println("n_level = $(n_level)")
=#
    return vert_level
end

bit_gdistances!(g::AbstractGraph{T}, source::Integer, vert_level::Vector{T}; queue_segment_size::Integer=20) where T<:Integer = 
bit_gdistances!(g, [source,], vert_level; queue_segment_size=20)

"""
    bit_gdistances(g, sources; queue_segment_size=20)
    bit_gdistances(g, source; queue_segment_size=20)

Parallel implementation of [`LightGraphs.gdistances!`](@ref) with dynamic load balancing.

### Optional Arguments
- `queue_segment_size = 20`: It is the number of vertices a thread can claim from a queue at a time.
For denser graphs, a smaller value of `queue_segment_size` could improve performance.

### References
- [Avoiding Locks and Atomic Instructions in Shared-Memory Parallel BFS Using Optimistic 
Parallelization](https://www.computer.org/csdl/proceedings/ipdpsw/2013/4979/00/4979b628-abs.html).
"""
bit_gdistances(g::AbstractGraph{T}, sources::Vector{<:Integer}; queue_segment_size::Integer=20) where T<:Integer = 
bit_gdistances!(g, sources, Vector{T}(undef, nv(g)); queue_segment_size=20)

bit_gdistances(g::AbstractGraph{T}, source::Integer; queue_segment_size::Integer=20) where T<:Integer = 
bit_gdistances!(g, [source,], Vector{T}(undef, nv(g)); queue_segment_size=20)

