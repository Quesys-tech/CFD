using Plots


G(x, β, z) = exp(-β * (x - z)^2)
F(x, α, a) = sqrt(max(1 - α^2 * (x - a)^2, 0))

function u_exact(x)
    z = -0.7
    δ = 0.005
    β = log(2) / 36δ^2
    a = 0.5
    α = 10
    if -0.8 <= x <= -0.6
        return 1 / 6 * (G(x, β, z - δ) + 4G(x, β, z) + G(x, β, z + δ))
    elseif -0.4 <= x <= -0.2
        return 1
    elseif 0 <= x <= 0.2
        return 1 - abs(10 * (x - 0.1))
    elseif 0.4 <= x <= 0.6
        return 1 / 6 * (F(x, α, a - δ) + 4F(x, α, a) + F(x, α, a + δ))
    else
        return 0
    end
end

function u_exact(x, t, c)
    y = x - c * t 
    z = ifelse(-1<=y<=1,y,mod(y+1,2)-1)
    u_exact(z)
end

function init(Δx::T) where {T}
    x = -1:Δx:1-Δx
    u = zeros(T, length(x))
    u .= u_exact.(x)
    x, u
end

function wenoz!(u⁺::Vector{T}, u::Vector{T}, c::T, Δx::T, Δt::T, ϵ=1e-10) where {T}
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

    @gif for n = 0:400
        plot(x, u, label="WENO")
        t = n * Δt
        plot!(x, u_exact.(x, t, c), label="Exact")
        wenoz!(u⁺, u, c, Δx, Δt)
        u .= u⁺
    end every 5
end

main()