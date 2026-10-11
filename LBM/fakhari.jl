module Fakhari

using LinearAlgebra
using StaticArrays
include("d2q9.jl")

export D2Q9Lattice, Parameters, PhaseFieldLBM, BGKCollision, MRTCollision,
    gamma, equilibria, relaxation_time, relaxation_matrix, apply_relaxation,
    local_stress, viscous_force, phase_derivatives!, macroscopic!, collide!,
    stream!, step!, diagnostics, interface_fronts, rti_parameters, initialize_rti

const CS2 = 1 / 3

"""
Fakhari parameters in lattice units; phi_L=0, phi_H=1.
tau is nu/cs2, NOT the conventional LBM relaxation time tau+1/2.
"""
struct Parameters{T <: AbstractFloat}
    rho_L::T
    rho_H::T
    xi::T
    sigma::T
    mobility::T
    tau_L::T
    tau_H::T
    gravity::SVector{2, T}
    viscosity_model::Symbol
    normal_tolerance::T
end

function Parameters(;
        T::Type{<:AbstractFloat} = Float64, rho_L = 1, rho_H = 3,
        xi = 5, sigma = 0, mobility = 0.02, tau_L = 0.1, tau_H = 0.1,
        gravity = (0, 0), viscosity_model = :linear,
        normal_tolerance = 1.0e-12
    )
    values = T.((rho_L, rho_H, xi, sigma, mobility, tau_L, tau_H, normal_tolerance))
    all(isfinite, values) || throw(ArgumentError("parameters must be finite"))
    rho_L > 0 && rho_H >= rho_L && xi > 0 && mobility > 0 &&
        tau_L > 0 && tau_H > 0 && sigma >= 0 && normal_tolerance >= 0 ||
        throw(ArgumentError("invalid density, interface, mobility or viscosity"))
    viscosity_model in (:linear, :harmonic, :dynamic) ||
        throw(ArgumentError("viscosity_model must be :linear, :harmonic or :dynamic"))
    g = SVector{2, T}(gravity)
    all(isfinite, g) || throw(ArgumentError("gravity must be finite"))
    return Parameters{T}(values[1:7]..., g, viscosity_model, values[8])
end

@inline function relaxation_time(p::Parameters{T}, phi, rho) where {T}
    if p.viscosity_model === :linear
        return p.tau_L + phi * (p.tau_H - p.tau_L)
    elseif p.viscosity_model === :harmonic
        return inv(inv(p.tau_L) + phi * (inv(p.tau_H) - inv(p.tau_L)))
    else
        cs2 = T(CS2)
        mu_L, mu_H = p.rho_L * cs2 * p.tau_L, p.rho_H * cs2 * p.tau_H
        return (mu_L + phi * (mu_H - mu_L)) / (rho * cs2)
    end
end

struct BGKCollision{T} <: AbstractCollision{2, 9, T} end

# d'Humieres / Lallemand-Luo orthogonal D2Q9 basis in OUR column ordering.
# Rows: density, energy, energy^2, jx, qx, jy, qy, pxx, pxy.
# Only the LAST TWO rows (trace-free normal stress and shear) have rate snu.
const MRT_BASIS = @SMatrix [
    1  1  1  1  1  1  1  1  1;
    -4 -1 -1 -1 -1  2  2  2  2;
    4 -2 -2 -2 -2  1  1  1  1;
    0  1  0 -1  0  1 -1 -1  1;
    0 -2  0  2  0  1 -1 -1  1;
    0  0  1  0 -1  1  1 -1 -1;
    0  0 -2  0  2  1  1 -1 -1;
    0  1 -1  1 -1  0  0  0  0;
    0  0  0  0  0  1 -1  1 -1
]

struct MRTCollision{T} <: AbstractCollision{2, 9, T}
    matrix::SMatrix{9, 9, T, 81}
    inverse::SMatrix{9, 9, T, 81}
end

function MRTCollision(::Type{T} = Float64) where {T}
    matrix = SMatrix{9, 9, T}(MRT_BASIS)
    return MRTCollision{T}(matrix, inv(matrix))
end

relaxation_matrix(::BGKCollision{T}, tau) where {T} =
    SMatrix{9, 9, T}(I) / (tau + T(0.5))
function relaxation_matrix(c::MRTCollision{T}, tau) where {T}
    s = inv(tau + T(0.5))
    rates = SVector{9, T}(1, 1, 1, 1, 1, 1, 1, s, s)
    return c.inverse * Diagonal(rates) * c.matrix
end

@inline apply_relaxation(::BGKCollision{T}, q, tau) where {T} =
    q / (tau + T(0.5))
@inline function apply_relaxation(c::MRTCollision{T}, q, tau) where {T}
    # K = I + (snu-1)*M^-1[:,8:9]*M[8:9,:], since all other rates are 1.
    # This is the full specified MRT operator, without a dense product per node.
    m8 = q[2] - q[3] + q[4] - q[5]
    m9 = q[6] - q[7] + q[8] - q[9]
    ds = inv(tau + T(0.5)) - one(T)
    return SVector{9, T}(
        ntuple(
            i ->
            q[i] + ds * (c.inverse[i, 8] * m8 + c.inverse[i, 9] * m9), 9
        )
    )
end

@inline function equilibria(lattice::D2Q9Lattice, phi, pstar, u)
    Gamma = gamma(lattice, u)
    return phi * Gamma, pstar * lattice.w + (Gamma - lattice.w)
end

"""Local second moment of K*(g-g_eq); never use g-gbar_eq here."""
@inline function local_stress(lattice::D2Q9Lattice{T}, collision, neq, tau) where {T}
    q = apply_relaxation(collision, neq, tau)
    xx = yy = xy = zero(T)
    @inbounds for i in 1:9
        ex, ey = lattice.c[1, i], lattice.c[2, i]
        xx += ex * ex * q[i]
        xy += ex * ey * q[i]
        yy += ey * ey * q[i]
    end
    return SMatrix{2, 2, T}(xx, xy, xy, yy)
end

@inline function viscous_force(lattice, collision, neq, tau, grad_rho)
    # -nu/(cs2*dt) * second_moment(K*neq) * grad(rho), nu=tau*cs2.
    # In BGK K=snu*I, giving -nu/((tau+1/2)*cs2*dt).
    return -tau * (local_stress(lattice, collision, neq, tau) * grad_rho)
end

mutable struct PhaseFieldLBM{T, C <: AbstractCollision{2, 9, T}}
    lattice::D2Q9Lattice{T}
    collision::C
    parameters::Parameters{T}
    h::Array{T, 3}             # (9, Nx, Ny), as in cavity.jl
    h_post::Array{T, 3}
    g::Array{T, 3}             # velocity-based distribution; sum(g)=pstar
    g_post::Array{T, 3}
    phi::Matrix{T}
    rho::Matrix{T}
    pstar::Matrix{T}
    u::Array{T, 3}             # (2, Nx, Ny)
    grad_phi::Array{T, 3}
    lap_phi::Matrix{T}
    tau::Matrix{T}
    force::Array{T, 3}
    force_mu::Array{T, 3}
    initial_mass::Float64
    iteration::Int
    threaded::Bool
end

# All kernels access only local g,u and neighboring phi. Rows can run independently.
function foreach_row!(f, s::PhaseFieldLBM)
    if s.threaded && Threads.nthreads() > 1
        Threads.@threads :static for y in axes(s.phi, 2)
            f(y)
        end
    else
        for y in axes(s.phi, 2)
            f(y)
        end
    end
    return nothing
end

function PhaseFieldLBM(
        phi::AbstractMatrix{T}, p::Parameters{T};
        u = zeros(T, 2, size(phi)...), pstar = zeros(T, size(phi)),
        collision = :mrt, threaded = true
    ) where {T}
    nx, ny = size(phi)
    nx >= 2 && ny >= 2 || throw(ArgumentError("domain needs at least 2x2 nodes"))
    size(u) == (2, nx, ny) && size(pstar) == (nx, ny) ||
        throw(DimensionMismatch("invalid u or pstar dimensions"))
    all(isfinite, phi) && all(isfinite, u) && all(isfinite, pstar) ||
        throw(ArgumentError("initial fields must be finite"))
    c = collision === :mrt ? MRTCollision(T) :
        collision === :bgk ? BGKCollision{T}() :
        throw(ArgumentError("collision must be :mrt or :bgk"))
    z2() = zeros(T, nx, ny)
    z3() = zeros(T, 2, nx, ny)
    z9() = zeros(T, 9, nx, ny)
    s = PhaseFieldLBM(
        D2Q9Lattice(T), c, p, z9(), z9(), z9(), z9(),
        Matrix{T}(phi), z2(), Matrix{T}(pstar), Array{T, 3}(u), z3(), z2(),
        z2(), z3(), z3(), sum(Float64, phi), 0, threaded
    )
    phase_derivatives!(s)
    material_fields!(s)
    all(>(0), s.rho) && all(>(0), s.tau) ||
        throw(ArgumentError("initial density and tau must be positive"))
    @inbounds for y in 1:ny, x in 1:nx
        heq, geq = equilibria(
            s.lattice, s.phi[x, y], s.pstar[x, y],
            SVector(s.u[1, x, y], s.u[2, x, y])
        )
        s.h[:, x, y] .= heq
        s.g[:, x, y] .= geq
    end
    forces_velocity!(s; recover_velocity = false)
    # Initialize transformed populations at gbar_eq and hbar_eq. The odd
    # half-force shift gives sum(g*e)=u-F/(2rho) and zero initial viscous stress.
    @inbounds for y in 1:ny, x in 1:nx
        fp = phase_forcing(s, x, y)
        fi = hydro_forcing(s, x, y)
        for i in 1:9
            s.h[i, x, y] -= fp[i] / 2
            s.g[i, x, y] -= fi[i] / 2
        end
    end
    return s
end

function phase_derivatives!(s::PhaseFieldLBM{T}) where {T}
    nx, ny = size(s.phi)
    inv_cs2 = T(3)
    foreach_row!(s) do y
        @inbounds for x in 1:nx
            gx = gy = lap = zero(T)
            center = s.phi[x, y]
            for i in 2:9
                xn = mod1(x + D2Q9_CX[i], nx)
                # Reflect phi at the halfway vertical walls (ghost=adjacent
                # fluid node): neutral, zero normal derivative. Together with
                # h bounce-back below this is our no-flux wall choice, NOT a
                # boundary formula attributed to Fakhari et al.
                yn = clamp(y + D2Q9_CY[i], 1, ny)
                neighbor = s.phi[xn, yn]
                weighted = s.lattice.w[i] * neighbor
                gx += D2Q9_CX[i] * weighted
                gy += D2Q9_CY[i] * weighted
                lap += s.lattice.w[i] * (neighbor - center)
            end
            s.grad_phi[1, x, y] = inv_cs2 * gx
            s.grad_phi[2, x, y] = inv_cs2 * gy
            s.lap_phi[x, y] = 2inv_cs2 * lap
        end
    end
    return s
end

function material_fields!(s::PhaseFieldLBM)
    p = s.parameters
    foreach_row!(s) do y
        @inbounds for x in axes(s.phi, 1)
            phi = s.phi[x, y]
            rho = p.rho_L + phi * (p.rho_H - p.rho_L)
            s.rho[x, y] = rho
            s.tau[x, y] = relaxation_time(p, phi, rho)
        end
    end
    return s
end

@inline function phase_forcing(s::PhaseFieldLBM{T}, x, y) where {T}
    p = s.parameters
    gx, gy = s.grad_phi[1, x, y], s.grad_phi[2, x, y]
    norm = hypot(gx, gy)
    normal = norm > p.normal_tolerance ? SVector(gx / norm, gy / norm) : zero(SVector{2, T})
    amplitude = (one(T) - 4(s.phi[x, y] - T(0.5))^2) / p.xi
    return SVector{9, T}(
        ntuple(
            i -> amplitude * s.lattice.w[i] *
                (s.lattice.c[1, i] * normal[1] + s.lattice.c[2, i] * normal[2]), 9
        )
    )
end

@inline function hydro_forcing(s::PhaseFieldLBM{T}, x, y) where {T}
    fx, fy = s.force[1, x, y], s.force[2, x, y]
    factor = inv(s.rho[x, y] * T(CS2))
    return SVector{9, T}(
        ntuple(
            i -> s.lattice.w[i] * factor *
                (s.lattice.c[1, i] * fx + s.lattice.c[2, i] * fy), 9
        )
    )
end

function forces_velocity!(s::PhaseFieldLBM{T}; recover_velocity = true) where {T}
    p = s.parameters
    beta, kappa = 12p.sigma / p.xi, T(1.5) * p.sigma * p.xi
    foreach_row!(s) do y
        @inbounds for x in axes(s.phi, 1)
            phi, rho, tau = s.phi[x, y], s.rho[x, y], s.tau[x, y]
            grad_phi = SVector(s.grad_phi[1, x, y], s.grad_phi[2, x, y])
            grad_rho = (p.rho_H - p.rho_L) * grad_phi
            chemical = 4beta * phi * (phi - one(T)) * (phi - T(0.5)) - kappa * s.lap_phi[x, y]
            Fs = chemical * grad_phi
            Fb = rho * p.gravity
            Fp = -s.pstar[x, y] * T(CS2) * grad_rho
            # Explicit time convention: streamed g^(n+1), phi^(n+1), pstar^(n+1),
            # but the regular equilibrium in the local stress uses u^n.
            # Then F^(n+1) recovers u^(n+1); collision uses this new velocity.
            # No pressure iteration or finite differences of velocity are used.
            u_lag = SVector(s.u[1, x, y], s.u[2, x, y])
            _, geq = equilibria(s.lattice, phi, s.pstar[x, y], u_lag)
            g = SVector{9, T}(ntuple(i -> s.g[i, x, y], 9))
            Fmu = viscous_force(s.lattice, s.collision, g - geq, tau, grad_rho)
            F = Fs + Fb + Fp + Fmu
            s.force_mu[1, x, y], s.force_mu[2, x, y] = Fmu
            s.force[1, x, y], s.force[2, x, y] = F
            if recover_velocity
                mx = my = zero(T)
                for i in 1:9
                    mx += D2Q9_CX[i] * g[i]
                    my += D2Q9_CY[i] * g[i]
                end
                s.u[1, x, y] = mx + F[1] / (2rho)
                s.u[2, x, y] = my + F[2] / (2rho)
            end
        end
    end
    return s
end

"""Reconstruct after streaming: phi=sum(h), pstar=sum(g), then F and u."""
function macroscopic!(s::PhaseFieldLBM{T}) where {T}
    foreach_row!(s) do y
        @inbounds for x in axes(s.phi, 1)
            phi = pressure = zero(T)
            for i in 1:9
                phi += s.h[i, x, y]
                pressure += s.g[i, x, y]
            end
            s.phi[x, y], s.pstar[x, y] = phi, pressure
        end
    end
    phase_derivatives!(s)
    material_fields!(s)
    forces_velocity!(s)
    return s
end

function collide!(s::PhaseFieldLBM{T}) where {T}
    tau_phi = s.parameters.mobility / T(CS2)
    omega_phi = inv(tau_phi + T(0.5))
    foreach_row!(s) do y
        @inbounds for x in axes(s.phi, 1)
            u = SVector(s.u[1, x, y], s.u[2, x, y])
            heq, geq = equilibria(s.lattice, s.phi[x, y], s.pstar[x, y], u)
            fp, fi = phase_forcing(s, x, y), hydro_forcing(s, x, y)
            g = SVector{9, T}(ntuple(i -> s.g[i, x, y], 9))
            omega = -apply_relaxation(s.collision, g - (geq - fi / 2), s.tau[x, y])
            for i in 1:9
                s.h_post[i, x, y] = s.h[i, x, y] -
                    omega_phi * (s.h[i, x, y] - (heq[i] - fp[i] / 2)) + fp[i]
                s.g_post[i, x, y] = g[i] + omega[i] + fi[i]
            end
        end
    end
    return s
end

function stream!(s::PhaseFieldLBM)
    nx, ny = size(s.phi)
    foreach_row!(s) do y
        @inbounds for x in 1:nx, i in 1:9
            xp = mod1(x - D2Q9_CX[i], nx)
            yp = y - D2Q9_CY[i]
            if 1 <= yp <= ny
                s.h[i, x, y] = s.h_post[i, xp, yp]
                s.g[i, x, y] = s.g_post[i, xp, yp]
            else
                # Stationary halfway bounce-back for g (no slip), and link
                # bounce-back for h (no phase crosses the wall). This maps
                # each outgoing population to one incoming population, so
                # streaming conserves phase mass, including diagonal links.
                opposite = D2Q9_OPPOSITE[i]
                s.h[i, x, y] = s.h_post[opposite, x, y]
                s.g[i, x, y] = s.g_post[opposite, x, y]
            end
        end
    end
    return s
end

function step!(s::PhaseFieldLBM, steps::Integer = 1)
    steps >= 0 || throw(ArgumentError("steps must be nonnegative"))
    for _ in 1:steps
        collide!(s)
        stream!(s)
        macroscopic!(s)
        s.iteration += 1
    end
    return s
end

function diagnostics(s::PhaseFieldLBM)
    mass = sum(Float64, s.phi)
    speed2 = maximum(
        s.u[1, x, y]^2 + s.u[2, x, y]^2
            for y in axes(s.phi, 2), x in axes(s.phi, 1)
    )
    finite = all(
        a -> all(isfinite, a), (
            s.h, s.g, s.phi, s.rho, s.pstar,
            s.u, s.grad_phi, s.lap_phi, s.tau, s.force, s.force_mu,
        )
    )
    return (
        step = s.iteration, phase_mass = mass, mass_drift = mass - s.initial_mass,
        relative_mass_drift = (mass - s.initial_mass) / max(abs(s.initial_mass), eps(Float64)),
        phi_min = minimum(s.phi), phi_max = maximum(s.phi),
        max_mach = sqrt(speed2 / CS2), finite = finite,
        rho_min = minimum(s.rho), tau_min = minimum(s.tau),
    )
end

function column_crossings(phi, x, y0)
    lo, hi, count = Inf, -Inf, 0
    @inbounds for y in 1:(size(phi, 2) - 1)
        a, b = phi[x, y] - 0.5, phi[x, y + 1] - 0.5
        if a == 0 || a * b < 0
            yi = y0 + (y - 1) + (a == 0 ? 0.0 : -a / (b - a))
            lo, hi, count = min(lo, yi), max(hi, yi), count + 1
        end
    end
    return count == 0 ? (NaN, NaN, 0) : (lo, hi, count)
end

"""phi=0.5 crossings: highest bubble at x=0, lowest heavy spike at x=L/2."""
function interface_fronts(s::PhaseFieldLBM; y0 = -size(s.phi, 2) / 2 + 0.5)
    nx = size(s.phi, 1)
    iseven(nx) || throw(ArgumentError("exact x=L/2 measurement needs even Nx"))
    _, bubble, nb = column_crossings(s.phi, 1, y0)
    spike, _, ns = column_crossings(s.phi, nx ÷ 2 + 1, y0)
    return (bubble = bubble, spike = spike, bubble_crossings = nb, spike_crossings = ns)
end

"""Fig.5 low-density-ratio RTI parameters, derived from dimensionless inputs."""
function rti_parameters(;
        T::Type{<:AbstractFloat} = Float64, L::Integer = 256,
        At = 0.5, Re = 3000, mu_star = 1, Ca = 0.26, Pe = 1000,
        t0 = 16000, rho_L = 1, xi = 5, viscosity_model = :linear
    )
    L > 0 && iseven(L) && 0 < At < 1 &&
        all(v -> isfinite(v) && v > 0, (Re, mu_star, Ca, Pe, t0, rho_L, xi)) ||
        throw(ArgumentError("invalid RTI dimensionless parameters"))
    rho_H = rho_L * (1 + At) / (1 - At)
    g = L / (At * t0^2)
    U0 = sqrt(g * L)
    mu_H = rho_H * U0 * L / Re
    mu_L = mu_H / mu_star
    nu_H, nu_L = mu_H / rho_H, mu_L / rho_L
    sigma, mobility = mu_H * U0 / Ca, U0 * L / Pe
    p = Parameters(;
        T, rho_L, rho_H, xi, sigma, mobility,
        tau_L = nu_L / CS2, tau_H = nu_H / CS2,
        gravity = (0, -g), viscosity_model
    )
    return (
        parameters = p, Nx = L, Ny = 4L, L = T(L), At = T(At), Re = T(Re),
        mu_star = T(mu_star), Ca = T(Ca), Pe = T(Pe), t0 = T(t0), g = T(g),
        U0 = T(U0), mu_H = T(mu_H), mu_L = T(mu_L), nu_H = T(nu_H), nu_L = T(nu_L),
        sigma = T(sigma), mobility = T(mobility), tau_phi = T(mobility / CS2),
    )
end

function initialize_rti(setup = rti_parameters(); collision = :mrt, threaded = true)
    p = setup.parameters
    T = typeof(p.rho_L)
    nx, ny, L = setup.Nx, setup.Ny, setup.L
    phi, pstar = zeros(T, nx, ny), zeros(T, nx, ny)
    @inbounds for x in 1:nx
        interface = T(0.1) * L * cos(2T(pi) * (x - 1) / L)
        for y in 1:ny
            yc = -2L + T(y) - T(0.5)
            # Initialization only: vertical d ≈ y-yI(x), NOT the exact
            # signed perpendicular distance of the paper's Eq.(40).
            d = yc - interface
            phi[x, y] = T(0.5) + T(0.5) * tanh(2d / p.xi)
        end
        # Columnwise hydrostatic physical pressure; p=0 at the top wall.
        # Convert to pstar=p/(rho*cs2). This is an initialization choice.
        pressure = zero(T)
        for y in ny:-1:1
            rho = p.rho_L + phi[x, y] * (p.rho_H - p.rho_L)
            pstar[x, y] = (pressure + setup.g * rho / 2) / (rho * T(CS2))
            pressure += setup.g * rho
        end
    end
    return PhaseFieldLBM(phi, p; pstar, collision, threaded)
end

end # module Fakhari
