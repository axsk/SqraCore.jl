## SQRA using VoronoiGraph.jl for the discretization

export sqra_voronoi

"""
    sqra_voronoi(u, beta, xs; nmc=1000)
Compute the voronoi diagram, approximate the volumes with `nmc` samples each and return the SQRA
"""
function sqra_voronoi(u, xs; beta=1, nmc=0)
    @assert length(u) == size(xs, 2)
    v, P = VoronoiGraph.voronoi(xs)
    if nmc > 0
        A, V = VoronoiGraph.mc_volumes(v, P, nmc)
    else
        A, V = VoronoiGraph.volumes(v, P)
    end
    C = sqra_weights(A, V, P)
    return sqra(u, C, beta=beta)
end

"""
    sqra_weights(A, Vs, P)

Given the boundary areas, volumes and centers, compute the SQRA weight matrix C
``C_ij = A_ij / (V_i * h_ij)``
"""
function sqra_weights(A, Vs, P)
    C = Diagonal(1 ./ Vs) * A
    divide_by_h!(C, P)

    # for unbounded boundary cells, set inbound rates to 0
    for (i, v) in enumerate(Vs)
        v < Inf && continue
        C[:, i] .= 0
    end

    return C
end

""" divide the areas `A[i,j]` by the the distance of its generators `P[i]-P[j]` """
function divide_by_h!(A, P)
    rows = rowvals(A)
    for j in axes(A, 2)
        for i in nzrange(A, j)
            i = rows[i]
            A[i,j] /= norm(P[i] - P[j])
        end
    end
end
