using StaticArrays
abstract type AbstractLBM end
abstract type AbstractD2LBM <: AbstractLBM end

struct D2Q9{T} <: AbstractD2LBM
    ω::SVector{9,T}# weight 9
    𝐜::SMatrix{2,9,T}
    ρ::Array{2,T} # density NxM
    𝐮::Array{3,T} # velocity 2xNxM
    ρ𝐮::Array{3,T} # momentum 2xNxM
    f⁰::Array{3,T} # distribution function before translation 9xNxM
    f::Array{3,T} # distribution function after translation 9xNxM

    function D2Q9(N, M, τ::T) where {T}
        ω = SA{T}[4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36]
        𝐜 = SA{T}[
            0 1 0 -1 0 1 -1 -1 1;
            0 0 1 0 -1 1 1 -1 -1
        ]
        ρ = Array{T}(undef, N, M)
        𝐮 = Array{T}(undef, 2, N, M)
        ρ𝐮 = Array{T}(undef, 2, N, M)
        f⁰ = Array{T}(undef, 9, N, M)
        f = Array{T}(undef, 9, N, M)
        new(ω, c, ρ, 𝐮)
    end
end

"""
    f_eq(ωₖ, 𝐜ₖ, 𝐮, ρ, c)

Compute the equilibrium distribution function.
"""
function f_eq(ωₖ, 𝐜ₖ, 𝐮, ρ, c)
    ωₖ * ρ * (1 + 3𝐜ₖ ⋅ 𝐮 / (c^2) + 9 * (𝐜ₖ ⋅ 𝐮)^2 / (2 * c^4) - 3 * 𝐮 ⋅ 𝐮 / (2 * c^2))
end

function density_velocity!(model::AbstractD2LBM)
    for i in axes(model.ρ𝐮, 2), j in axes(model.ρ𝐮, 3)
        model.ρ[i, j] = 0
        model.ρ𝐮[:, i, j] .= 0
        for (k, fₖ) in enumerate(model.f[:, i, j])
            model.ρ[i, j] += fₖ
            @. model.ρ𝐮[:, i, j] += fₖ * model.𝐜[:, k]
        end
    end
    for i in axes(model.ρ𝐮, 2), j in axes(model.ρ𝐮, 3)
        model.𝐮[:, i, j] .= model.ρ𝐮[:, i, j] ./ model.ρ[i, j]
    end
end