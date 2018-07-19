@testset "Generate Map" begin

    function make_vec(g::AbstractGraph{T}) where T<:Integer
        return Vector{T}(undef, nv(g))
    end

    g1 = StarGraph(5)

    for g in testgraphs(g1)
        s = @inferred(LightGraphs.generate_min_colors(g, greedy_color, 5))
        @test s.num_colors == 2
    end

    for g in testgraphs(g1)
        s = @inferred(LightGraphs.generate_max_set(g, make_vec, 5))
        @test length(s) == 5
    end

    for g in testgraphs(g1)
        s = @inferred(LightGraphs.generate_max_set(g, make_vec, 5))
        @test length(s) == 5
    end
end
