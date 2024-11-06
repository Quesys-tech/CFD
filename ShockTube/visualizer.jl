using Plots
using ForwardDiff
using CSV
using DataFrames

function exact(ρₗ, pₗ, uₗ, ρᵣ, pᵣ, uᵣ, t, x; γ=1.4)

    ρ = zeros(length(x))
    p = zeros(length(x))
    u = zeros(length(x))

    if t <= 0
        for i in eachindex(x)
            ρ[i] = ifelse(x[i] < 0, ρₗ, ρᵣ)
            p[i] = ifelse(x[i] < 0, pₗ, pᵣ)
            u[i] = ifelse(x[i] < 0, uₗ, uᵣ)
        end
    else
        # solve pressure ratio
        cᵣ = sqrt(γ * pᵣ / ρᵣ)
        cₗ = sqrt(γ * pₗ / ρₗ)

        function pratioerr(P)
            lhs = sqrt(2 / γ) * (P - 1) / sqrt(γ - 1 + (γ + 1) * P)
            rhs = 2 / (γ - 1) * cₗ / cᵣ * (1 - (pᵣ / pₗ * P)^((γ - 1) / 2γ))
            rhs += (uₗ - uᵣ) / cᵣ
            return lhs - rhs
        end

        P = pₗ
        for i = 1:10
            abs(pratioerr(P)) < 1e-15 && break
            P -= pratioerr(P) / ForwardDiff.derivative(pratioerr, P)
        end

        p₂ = pᵣ * P
        ρ₂ = ρᵣ * (P + (γ - 1) / (γ + 1)) / ((γ - 1) / (γ + 1) * P + 1)
        u₂ = uᵣ + cᵣ * sqrt(2 / γ) * (P - 1) / sqrt(γ - 1 + (γ + 1) * P)
        V_s = uᵣ + (P - 1) * cᵣ^2 / (γ * (u₂ - uᵣ))

        V_c = u₂
        V_rt = u₂ - ((γ - 1) / 2 * (uₗ - u₂) + cₗ)
        V_rh = uₗ - cₗ

        for i in eachindex(x)
            if x[i] > V_s * t # front of shock wave
                ρ[i] = ρᵣ
                p[i] = pᵣ
                u[i] = uᵣ
            elseif x[i] > V_c * t # shock 
                ρ[i] = ρ₂
                p[i] = p₂
                u[i] = u₂
            elseif x[i] > V_rt * t # contact discontinuity
                ρ₃ = ρₗ * (p₂ / pₗ)^(1 / γ)
                ρ[i] = ρ₃
                p[i] = p₂
                u[i] = u₂
            elseif x[i] > V_rh * t # expansion fan
                u₄ = 2 / (γ + 1) * (x[i] / t + cₗ + (γ - 1) / 2 * uₗ)
                c₄ = cₗ - (γ - 1) / 2 * (u₄ - uₗ)
                p₄ = pₗ * (c₄ / cₗ)^(2 * γ / (γ - 1))
                ρ₄ = ρₗ * (p₄ / pₗ)^(1 / γ)
                ρ[i] = ρ₄
                p[i] = p₄
                u[i] = u₄
            else
                ρ[i] = ρₗ
                p[i] = pₗ
                u[i] = uₗ
            end
        end
    end
    return ρ, p, u
end

function main()
    x = -1:0.005:1
    ρ_exact, p_exact, u_exact = exact(1, 1, 0, 0.125, 0.1, 0, 0.400208016927473, x)
    df = CSV.read("sod.csv", DataFrame)

    # plot density
    plot(x, ρ_exact, label="exact", xlabel="x", ylabel="ρ", legend=:topleft)
    plot!(df.x, df.density, label="numerical")
    savefig("density.png")

    # plot velocity
    plot(x, u_exact, label="exact", xlabel="x", ylabel="u", legend=:topleft)
    plot!(df.x, df.velocity, label="numerical")
    savefig("velocity.png")

    # plot pressure
    plot(x, p_exact, label="exact", xlabel="x", ylabel="p", legend=:topleft)
    plot!(df.x, df.pressure, label="numerical")
    savefig("pressure.png")

end
main()