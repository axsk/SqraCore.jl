using SqraCore

# using a regular grid
sqra_grid(rand(10))
sqra_grid(rand(5,5,5), beta=3, h=1/5)

# using voronoi tesselations
using VoronoiGraph
U(x::Vector) = sum(x.^2)
U(x) = mapslices(U, x, dims=1) |> vec
xs = rand(3,100)

sqra_voronoi(U(xs), xs)

# TODO: test if the results are correct
