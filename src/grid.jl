### SQRA on regular grids

export sqra_grid

""" Convenience wrapper for the SQRA,
implicitly computing the connectivity for the N-D Array `u` based on a regular grid """
function sqra_grid(u::Array; beta::Real=1, h=1)
    A = grid_adjacency(size(u); h)
    sqra(vec(u), A, beta=beta)
end

function sqra_grid(u::Vector; kwargs...)
    u = reshape(u, :, 1)
    sqra_grid(u; kwargs...)
end

""" compute the adjacency matrix for a regular grid
in multiple dimensions with uniform cell-length h"""
function grid_adjacency(dims::NTuple{N, Int} where N; h=1)
    dims = reverse(dims) # somehow we have to take the krons backwards, dont know why
    k = []
    for i in eachindex(dims)
        x = [spdiagm(ones(s)) for s in dims] # identity in all dimensions
        x[i] = spdiagm(-1 => ones(dims[i]-1), 1=>ones(dims[i]-1)) # neighbour matrix in dimension i
        push!(k, kron(x...))
    end
    sum(k) ./ h^2
end
