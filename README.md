Implementation of approximation algorithms for the travelling salesman problem using LightGraphs.jl library.

src/travellingsalesman/metric_travelling_salesman.jl is an approximation algorithm for the case where edge weights follow the triangle inequality. The algorithm can be found in sub-section 2.1 of the following [post](http://www.cs.tufts.edu/~cowen/advanced/2002/adv-lect3.pdf).

Corresponding tests were added in `test/travellingsalesman/metric_travelling_salesman.jl`.
