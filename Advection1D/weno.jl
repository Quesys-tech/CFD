using Plots

function u_exact(x, t, c)
    y = x - c * t
    while !(0 <= y <= 1)
        if y < 0
            y += 1
        else
            y -= 1
        end
    end
    if 1/3 <= y <= 2/3
        return 1
    else
        return 0
    end
end

function init(Δx::T) where {T}
    x = 0:Δx:1-Δx
    u = zeros(T, length(x))
    u .= u_exact.(x, 0, 0)
    x, u
end

function weno_pbc!(u⁺::Vector{T}, u::Vector{T}, c::T, Δx::T, Δt::T, ϵ=1e-10) where T
    N = length(u)
    @assert length(u⁺) == N
    @assert c > 0
    @assert Δx > 0
    @assert Δt > 0

    function D⁻(i)
        uᵢ = u[mod1(i, N)]
        uᵢ₋₁ = u[mod1(i - 1, N)]
        (uᵢ - uᵢ₋₁) / Δx
    end

    for i in eachindex(u)
        D⁻ᵢ₋₂ = D⁻(i - 2)
        D⁻ᵢ₋₁ = D⁻(i - 1)
        D⁻ᵢ = D⁻(i)
        D⁻ᵢ₊₁ = D⁻(i + 1)
        D⁻ᵢ₊₂ = D⁻(i + 2)

        b₁ = D⁻ᵢ₋₂ / 3 - 7D⁻ᵢ₋₁ / 6 + 11D⁻ᵢ / 6
        b₂ = -D⁻ᵢ₋₁ / 6 + 5D⁻ᵢ / 6 + D⁻ᵢ₊₁ / 3
        b₃ = D⁻ᵢ / 3 + 5D⁻ᵢ₊₁ / 6 - D⁻ᵢ₊₂ / 6
        S₁ = 13 / 12 * (D⁻ᵢ₋₂ - 2D⁻ᵢ₋₁ + D⁻ᵢ)^2 + 1 / 4 * (D⁻ᵢ₋₂ - 4D⁻ᵢ₋₁ + 3D⁻ᵢ)^2
        S₂ = 13 / 12 * (D⁻ᵢ₋₁ - 2D⁻ᵢ + D⁻ᵢ₊₁)^2 + 1 / 4 * (D⁻ᵢ₋₁ - D⁻ᵢ)^2
        S₃ = 13 / 12 * (D⁻ᵢ - 2D⁻ᵢ₊₁ + D⁻ᵢ₊₂)^2 + 1 / 4 * (3D⁻ᵢ - 4D⁻ᵢ₊₁ + D⁻ᵢ₊₂)^2
        α₁ = 0.1 / (S₁ + ϵ)^2
        α₂ = 0.6 / (S₂ + ϵ)^2
        α₃ = 0.3 / (S₃ + ϵ)^2
        sum_α = α₁ + α₂ + α₃
        w₁ = α₁ / sum_α
        w₂ = α₂ / sum_α
        w₃ = α₃ / sum_α
        u⁺[i] = u[i] - Δt * c * (w₁ * b₁ + w₂ * b₂ + w₃ * b₃)
    end
end

function main()
    c = 1.0
    Δx = 0.01
    Δt = 0.5Δx / c

    ν = c * Δt / Δx
    @show ν

    x, u = init(Δx)
    u⁺ = similar(u)

    @gif for n = 0:200
        plot(x, u, label="WENO")
        t = n * Δt
        plot!(x, u_exact.(x, t, c), label="Exact",title="\$t=$(round(t,digits=4))\$")
        weno_pbc!(u⁺, u, c, Δx, Δt)
        u .= u⁺
    end
end

main()