export bit_gdistances, bit_gdistances!

function bit_gdistances!(g::AbstractGraph{T}, source, vert_level; sort_alg = QuickSort) where T
    n = nv(g)
    visited = falses(n)
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

#count_E = count_bit = count_q = 0
    while !isempty(cur_level)
        @inbounds for v in cur_level
#count_q+=1
            @inbounds @simd for i in outneighbors(g, v)
#count_E+=1
                if !visited[i]
#count_bit+=1
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
#println("count_bit = count_bit")
    return vert_level
end

bit_gdistances(g::AbstractGraph{T}, source; sort_alg = Base.Sort.QuickSort) where T = bit_gdistances!(g, source, fill(typemax(T), nv(g)); sort_alg = sort_alg)

