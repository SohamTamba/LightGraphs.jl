@testset "Johnson" begin
    g3 = PathGraph(5)
    d = [0 1 2 3 4; 5 0 6 7 8; 9 10 0 11 12; 13 14 15 0 16; 17 18 19 20 0]
    for g in testgraphs(g3)
        z3 = @inferred(johnson_shortest_paths(g, d, parallel = true))
        @test z3.dists[1, :] == [0, 1, 7, 18, 34]
        @test z3.parents[3, :] == [2, 3, 0, 3, 4]

        @test @inferred(enumerate_paths(z3))[2][2] == []
        @test @inferred(enumerate_paths(z3))[2][4] == [2, 3, 4]
    end

    g4 = PathDiGraph(4)
    d = -ones(4, 4)
    for g in testdigraphs(g4)
        z4 = @inferred(johnson_shortest_paths(g, d))
        @test length(enumerate_paths(z4, 4, 3)) == 0
        @test length(enumerate_paths(z4, 4, 1)) == 0
        @test length(enumerate_paths(z4, 2, 3)) == 2
        @test z4.dists[1, :] == [0, -1, -2, -3]
        @test z4.parents[3, :] == [0, 0, 0, 3]
    end

    g5 = PathGraph(5)
    d = -ones(5, 5)
    for g in testgraphs(g5)
        @test_throws LightGraphs.NegativeCycleError johnson_shortest_paths(g, d)
    end
end