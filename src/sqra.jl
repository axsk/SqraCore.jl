using SparseArrays
using LinearAlgebra

export sqra

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
