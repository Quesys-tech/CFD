using StaticArrays

abstract type AbstractLattice{D, Q, T} end
abstract type AbstractCollision{D, Q, T} end

# Direction-first layout: (rest, E, N, W, S, NE, NW, SW, SE).
const D2Q9_CX = (0, 1, 0, -1, 0, 1, -1, -1, 1)
const D2Q9_CY = (0, 0, 1, 0, -1, 1, 1, -1, -1)
const D2Q9_OPPOSITE = (1, 4, 5, 2, 3, 8, 9, 6, 7)

struct D2Q9Lattice{T} <: AbstractLattice{2, 9, T}
    w::SVector{9, T}
    c::SMatrix{2, 9, T, 18}
    function D2Q9Lattice(::Type{T}) where {T}
        w = @SVector [4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36]
        c = @SMatrix [
            0 1 0 -1 0 1 -1 -1 1;
            0 0 1 0 -1 1 1 -1 -1
        ]
        return new{T}(w, c)
    end
end

"""Unit-density second-order equilibrium, with dx = dt = c = 1."""
@inline function gamma_i(w::T, c, u) where {T}
    eu = c[1] * u[1] + c[2] * u[2]
    u2 = u[1]^2 + u[2]^2
    return w * (one(T) + 3eu + T(4.5) * eu^2 - T(1.5) * u2)
end

@inline gamma(lattice::D2Q9Lattice, u) =
    SVector(ntuple(i -> gamma_i(lattice.w[i], lattice.c[:, i], u), 9))
