using SparseArrays
using LinearAlgebra


### SQRA core routine

""" Compute the Square-Root approximation to the generator of Langevin diffusion
for a given potential vector `u` and connectivity matrix `C` at coldness `beta`.
`C` describes the geometry via C_ij = A_ij / (h_ij * V_i)
with A, V, and h the area, volume and distances.
"""
function sqra(u::Vector, C::SparseMatrixCSC; beta::Real)
    # We use sqrt(p_j/p_i) = exp(-beta/2 (U_j - U_i)) where the latter is num. more stable
    I, J, a = findnz(C)
    v = -beta / 2 .* u  # TODO: do we need this intermediate allocation?
    q = similar(u, size(a))
    for n in eachindex(q)
        q[n] = exp(v[J[n]] - v[I[n]]) * a[n] / beta
    end
    Q = sparse(I, J, q, size(C)...)
    Q = setdiagonal(Q)
end

function setdiagonal(Q)
	Q = Q - spdiagm(diag(Q)) # remove diagonal
    Q = Q - spdiagm(sum(Q, dims=2)|>vec) # rowsum 0
	return Q
end


### SQRA on regular grids

""" Convenience wrapper for the SQRA,
implicitly computing the connectivity for the N-D Array `u` based on a regular grid """
function sqra(u::Array; beta::Real=1)
    A = grid_adjacency(size(u))
    sqra(vec(u), A, beta=beta)
end

function sqra(u::Vector; beta::Real=1)
    u = reshape(u, :, 1)
    sqra(u, beta=beta)
end

""" compute the adjacency matrix for a regular grid in multiple dimensions with Δ=1 """
function grid_adjacency(dims::NTuple{N, Int} where N)
    dims = reverse(dims) # somehow we have to take the krons backwards, dont know why
    k = []
    for i in eachindex(dims)
        x = [spdiagm(ones(s)) for s in dims] # identity in all dimensions
        x[i] = spdiagm(-1 => ones(dims[i]-1), 1=>ones(dims[i]-1)) # neighbour matrix in dimension i
        push!(k, kron(x...))
    end
    sum(k)
end


### SQRA on voronoi cells

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
