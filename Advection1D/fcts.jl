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
    if 0.25 <= y <= 0.75
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

function ftcs_pbc!(u⁺::Vector{T}, u::Vector{T}, ν) where {T}
    @assert size(u⁺) == size(u)

    for i in eachindex(u)
        uᵣ::T = 0
        uₗ::T = 0
        if i == 1
            uₗ = u[lastindex(u)]
            uᵣ = u[i+1]
        elseif i == lastindex(u)
            uₗ = u[i-1]
            uᵣ = u[1]
        else
            uₗ = u[i-1]
            uᵣ = u[i+1]
        end
        u⁺[i] = u[i] - 0.5ν * (uᵣ - uₗ)
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

    @gif for n=0:200
        plot(x,u,label="FCTS")
        t = n*Δt
        plot(x,u_exact.(x,t,c),label="Exact")
        ftcs_pbc!(u⁺,u,ν)
        u.= u⁺
    end
end 

main()