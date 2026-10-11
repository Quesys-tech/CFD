using DelimitedFiles
using DataFrames
using CSV

const γ = 1.4
const CFL = 0.2
const Δx = 0.005

function setup_sod(x)
    Q = zeros(3, length(x))

    for i in eachindex(x)
        if x[i] < 0
            ρ = 1.0
            u = 0.0
            p = 1.0
            Q[1, i] = ρ
            Q[2, i] = ρ * u
            Q[3, i] = p / (γ - 1) + ρ * u^2 / 2
        else
            ρ = 1 / 8
            u = 0.0
            p = 0.1
            Q[1, i] = ρ
            Q[2, i] = ρ * u
            Q[3, i] = p / (γ - 1) + ρ * u^2 / 2
        end
    end
    return Q
end

p(Q) = (γ - 1) * (Q[3] - Q[2]^2 / 2Q[1])
a(Q) = sqrt(γ * p(Q) / Q[1])
# H(Q) = (Q[3] + p(Q)) / Q[1]

function flux_ausm!(𝐔ⱼ, 𝐔ⱼ₊₁)
    M⁺(M) = ifelse(abs(M) > 1, (M + abs(M)) / 2, (M + 1)^2 / 4)
    M⁻(M) = ifelse(abs(M) > 1, (M - abs(M)) / 2, -(M - 1)^2 / 4)
    Mⱼ = 𝐔ⱼ[2] / 𝐔ⱼ[1] / a(𝐔ⱼ)
    Mⱼ₊₁ = 𝐔ⱼ₊₁[2] / 𝐔ⱼ₊₁[1] / a(𝐔ⱼ₊₁)
    m₁₂ = M⁺(Mⱼ) + M⁻(Mⱼ₊₁)

    function 𝚽(𝐔)
        return [𝐔[1], 𝐔[2], 𝐔[3] + p(𝐔)] .* a(𝐔)
    end
    𝚽₁₂ = ifelse(m₁₂ ≥ 0, 𝚽(𝐔ⱼ), 𝚽(𝐔ⱼ₊₁))
    𝐟ᶜ₁₂ = m₁₂ .* 𝚽₁₂

    P⁺(M) = ifelse(abs(M) > 1, (1 + sign(M)) / 2, (M + 1)^2 * (2 - M) / 4)
    P⁻(M) = ifelse(abs(M) > 1, (1 - sign(M)) / 2, (M - 1)^2 * (2 + M) / 4)
    p₁₂ = P⁺(Mⱼ) * p(𝐔ⱼ) + P⁻(Mⱼ₊₁) * p(𝐔ⱼ₊₁)
    𝐩₁₂ = [0, p₁₂, 0]
    𝐟₁₂ = 𝐟ᶜ₁₂ .+ 𝐩₁₂
    return 𝐟₁₂
end

function sod_solver(x)
    𝐔 = setup_sod(x)
    𝐅 = zeros(size(𝐔, 1), size(𝐔, 2) - 1)
    t = 0

    while t < 0.4
        maxu = maximum(𝐔[2, :] ./ 𝐔[1, :] .+ sqrt.(γ * p(𝐔) ./ 𝐔[1, :]))
        Δt = CFL * Δx / maxu

        for i in 1:(length(x) - 1)
            𝐅[:, i] = flux_ausm!(𝐔[:, i], 𝐔[:, i + 1])
        end
        for i in 2:(length(x) - 1)
            𝐔[:, i] .-= (Δt / Δx) .* (𝐅[:, i] .- 𝐅[:, i - 1])
        end
        t += Δt
    end
    @info "t = $t"
    return 𝐔
end

function main()
    x = -1:Δx:1
    @time 𝐔 = sod_solver(x)
    @time 𝐔 = sod_solver(x)

    ρ = 𝐔[1, :]
    u = 𝐔[2, :] ./ 𝐔[1, :]

    df = DataFrame(x = x, density = ρ, velocity = u, pressure = [p(𝐔[:, i]) for i in eachindex(x)])
    return CSV.write("sod.csv", df)
end

main()
