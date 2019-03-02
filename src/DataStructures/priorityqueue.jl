# This file contains code that was formerly a part of Julia. License is MIT: http://julialang.org/license

# Int_PriorityQueue
# -------------

module IPQ

import DataStructures: peek, dequeue!, dequeue_pair!
import Base: setindex!, getindex, length, lt, isempty
using DataStructures: Ordering, ForwardOrdering, ReverseOrdering, Forward, Reverse


export Int_PriorityQueue

heapparent(i::K) where {K<:Integer} = div(i, 2)
heapleft(i::K) where {K<:Integer} = 2*i
heapright(i::K) where {K<:Integer} = 2*i+1





mutable struct Int_PriorityQueue{K<:Integer, V, O<:Ordering}
    # Binary heap of (element, priority) pairs.
    xs::Array{Pair{K,V}, 1}
    o::O

    # Map elements to their index in xs
    index::Array{K, 1}    
    n::K

    function Int_PriorityQueue{K,V}(n::Int) where K<:Integer where V
        new{K,V,Base.Order.ForwardOrdering}(Vector{Pair{K,V}}(), Base.Order.Forward, zeros(K, n), K(n))
    end
end




length(pq::Int_PriorityQueue) = length(pq.xs)
isempty(pq::Int_PriorityQueue) = isempty(pq.xs)
haskey(pq::Int_PriorityQueue, key) = (key <=pq.n && key > 0 && pq.index[key] != 0)


peek(pq::Int_PriorityQueue) = pq.xs[1]

function percolate_down!(pq::Int_PriorityQueue, i::Integer)
    x = pq.xs[i]
    @inbounds while (l = heapleft(i)) <= length(pq)
        r = heapright(i)
        j = r > length(pq) || lt(pq.o, pq.xs[l].second, pq.xs[r].second) ? l : r
        if lt(pq.o, pq.xs[j].second, x.second)
            pq.index[pq.xs[j].first] = i
            pq.xs[i] = pq.xs[j]
            i = j
        else
            break
        end
    end
    pq.index[x.first] = i
    pq.xs[i] = x
end


function percolate_up!(pq::Int_PriorityQueue, i::Integer)
    x = pq.xs[i]
    @inbounds while i > 1
        j = heapparent(i)
        if lt(pq.o, x.second, pq.xs[j].second)
            pq.index[pq.xs[j].first] = i
            pq.xs[i] = pq.xs[j]
            i = j
        else
            break
        end
    end
    pq.index[x.first] = i
    pq.xs[i] = x
end

# Equivalent to percolate_up! with an element having lower priority than any other
function force_up!(pq::Int_PriorityQueue, i::Integer)
    x = pq.xs[i]
    @inbounds while i > 1
        j = heapparent(i)
        pq.index[pq.xs[j].first] = i
        pq.xs[i] = pq.xs[j]
        i = j
    end
    pq.index[x.first] = i
    pq.xs[i] = x
end

function getindex(pq::Int_PriorityQueue{K,V}, key) where {K<:Integer,V}
    pq.xs[pq.index[key]].second
end


function get(pq::Int_PriorityQueue{K,V}, key, deflt) where {K<:Integer,V}
    i = pq.index[key]
    i == 0 ? deflt : pq.xs[i].second
end


# Change the priority of an existing element, or equeue it if it isn't present.
function setindex!(pq::Int_PriorityQueue{K, V}, value, key) where {K<:Integer,V}
    if !haskey(pq, key)
        if key > pq.n
            resize!(pq.index, key)
            pq.index[pq.n+1:key] .= zero(K)
            pq.n = key
        end
        push!(pq.xs, Pair{K, V}(key, value))
        pq.index[key] = length(pq.xs)
    end

    i = pq.index[key]
    oldvalue = pq.xs[i].second
    pq.xs[i] = Pair{K,V}(key, value)
    if lt(pq.o, oldvalue, value)
        percolate_down!(pq, i)
    else
        percolate_up!(pq, i)
    end
    value
end


function dequeue!(pq::Int_PriorityQueue{K, V}) where K<:Integer where V<:Integer
    x = pq.xs[1]
    y = pop!(pq.xs)
    if !isempty(pq)
        pq.xs[1] = y
        pq.index[y.first] = 1
        percolate_down!(pq, 1)
    end
    pq.index[x.first] = zero(K)
    x.first
end

function dequeue!(pq::Int_PriorityQueue, key)
    idx = pq.index[key]
    force_up!(pq, idx)
    dequeue!(pq)
    key
end

function dequeue_pair!(pq::Int_PriorityQueue)
    x = pq.xs[1]
    y = pop!(pq.xs)
    if !isempty(pq)
        pq.xs[1] = y
        pq.index[y.first] = 1
        percolate_down!(pq, 1)
    end
    delete!(pq.index, x.first)
    x
end

function dequeue_pair!(pq::Int_PriorityQueue, key)
    idx = pq.index[key]
    force_up!(pq, idx)
    dequeue_pair!(pq)
end

function delete!(pq::Int_PriorityQueue, key)
    dequeue_pair!(pq, key)
    pq
end

end
