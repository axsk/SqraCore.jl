module SqraCore

include("sqra.jl")
include("grid.jl")
include("utils.jl")

# optional code binding to VoronoiGraph
using Requires
function __init__()
    @require VoronoiGraph="e5f3e3e8-00d7-4f74-8011-648b521326aa" include("voronoi.jl")
end

end
