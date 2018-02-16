using BenchmarkTools, Compat, LightGraphs, SimpleWeightedGraphs

srand(7)

m = 4
g = 0

for n in [1000, 10000, 100000]

	g = SimpleWeightedGraph(n)

	for i in 2:n
		add_edge!(g, i-1, i, 2*n)
	end 

	for i in 1:m*n+1
		add_edge!(g, rand(1:n), rand(1:n), rand(1:2*n))
	end
    
    print("For |V| = $n, \nOld Implementation: ")
    @btime kruskal_mst_old(g)
	print("Optimized Implementation: ")
	@btime kruskal_mst(g)
end