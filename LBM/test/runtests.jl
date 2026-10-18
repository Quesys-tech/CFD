using Test
using LinearAlgebra
using StaticArrays
using Random
include("../fakhari.jl")
using .Fakhari

# Independent dense definitions, used to check the optimized collision and stress.
function dense_K(collision, tau)
    s = 1 / (tau + 0.5)
    collision isa BGKCollision && return Matrix{Float64}(I, 9, 9) * s
    M = Matrix(collision.matrix)
    return inv(M) * Diagonal([1, 1, 1, 1, 1, 1, 1, s, s]) * M
end

function direct_equilibria(lattice, phi, pressure, u)
    Gamma = [
        lattice.w[i] * (
            1 + dot(lattice.c[:, i], u) / Fakhari.CS2 +
                dot(lattice.c[:, i], u)^2 / (2Fakhari.CS2^2) - dot(u, u) / (2Fakhari.CS2)
        ) for i in 1:9
    ]
    return phi * Gamma, pressure * lattice.w + Gamma - lattice.w
end

@testset "D2Q9 ordering and equilibrium moments" begin
    for T in (Float64, Float32)
        lattice = D2Q9Lattice(T)
        tol = 10eps(T)
        @test sum(lattice.w) ≈ 1 atol = tol
        @test lattice.c * lattice.w ≈ zeros(T, 2) atol = tol
        @test lattice.c * Diagonal(lattice.w) * lattice.c' ≈ Fakhari.CS2 * I atol = tol
        for i in 1:9
            @test lattice.c[:, Fakhari.D2Q9_OPPOSITE[i]] == -lattice.c[:, i]
        end
        u = SVector{2, T}(0.034, -0.017)
        phi, pressure = T(0.37), T(0.82)
        Gamma = gamma(lattice, u)
        heq, geq = equilibria(lattice, phi, pressure, u)
        @test sum(Gamma) ≈ 1 atol = tol
        @test lattice.c * Gamma ≈ u atol = tol
        @test sum(heq) ≈ phi atol = tol
        @test sum(geq) ≈ pressure atol = tol
        @test lattice.c * geq ≈ u atol = tol
        @test eltype(heq) == T
        @test eltype(geq) == T
        reference_h, reference_g = direct_equilibria(lattice, phi, pressure, u)
        @test heq ≈ reference_h atol = tol
        @test geq ≈ reference_g atol = tol
    end
end

@testset "Isotropic phi stencil, including walls and periodic seam" begin
    p = Parameters()
    s = PhaseFieldLBM(fill(0.37, 16, 32), p; threaded = false)
    @test maximum(abs, s.grad_phi) < 1.0e-15
    @test maximum(abs, s.lap_phi) < 1.0e-15
    # An affine field has exact centered derivatives away from the walls/seam.
    s.phi .= [0.1 + 0.001x + 0.002y for x in 1:16, y in 1:32]
    phase_derivatives!(s)
    @test s.grad_phi[1, 2:15, 2:31] ≈ fill(0.001, 14, 30) atol = 1.0e-15
    @test s.grad_phi[2, 2:15, 2:31] ≈ fill(0.002, 14, 30) atol = 1.0e-15
    @test maximum(abs, s.lap_phi[2:15, 2:31]) < 1.0e-15
    # Independently evaluate the specified lattice sums on ALL nodes.
    nx, ny = size(s.phi)
    for y in 1:ny, x in 1:nx
        grad = zeros(2)
        lap = 0.0
        for i in 1:9
            ex, ey = Int.(s.lattice.c[:, i])
            neighbor = s.phi[mod1(x + ex, nx), clamp(y + ey, 1, ny)]
            grad += s.lattice.w[i] * s.lattice.c[:, i] * neighbor / Fakhari.CS2
            lap += 2s.lattice.w[i] * (neighbor - s.phi[x, y]) / Fakhari.CS2
        end
        @test s.grad_phi[:, x, y] ≈ grad atol = 1.0e-15
        @test s.lap_phi[x, y] ≈ lap atol = 1.0e-15
    end
end

@testset "MRT modes and local non-equilibrium stress (BGK and MRT)" begin
    lattice = D2Q9Lattice(Float64)
    rng = MersenneTwister(17)
    q = SVector{9}(randn(rng, 9))
    for collision in (BGKCollision{Float64}(), MRTCollision())
        for tau in (0.0058, 0.0174, 0.4)
            K = dense_K(collision, tau)
            @test Matrix(relaxation_matrix(collision, tau)) ≈ K atol = 1.0e-14
            @test apply_relaxation(collision, q, tau) ≈ K * q atol = 1.0e-14
            reference_stress = zeros(2, 2)
            for i in 1:9
                reference_stress += (lattice.c[:, i] * lattice.c[:, i]') * (K * q)[i]
            end
            grad_rho = SVector(0.021, -0.043)
            @test local_stress(lattice, collision, q, tau) ≈ reference_stress atol = 1.0e-14
            @test viscous_force(lattice, collision, q, tau, grad_rho) ≈
                -tau * reference_stress * grad_rho atol = 1.0e-14
            if collision isa MRTCollision
                M = Matrix(collision.matrix)
                @test M * Matrix(relaxation_matrix(collision, tau)) * inv(M) ≈
                    Diagonal([ones(7); 1 / (tau + 0.5); 1 / (tau + 0.5)]) atol = 1.0e-14
                @test M[8, :] == vec(lattice.c[1, :] .^ 2 - lattice.c[2, :] .^ 2)
                @test M[9, :] == vec(lattice.c[1, :] .* lattice.c[2, :])
            end
        end
    end
end

@testset "Forces, independent pressure/velocity reconstruction and collision" begin
    nx, ny = 10, 20
    p = Parameters(
        rho_L = 1, rho_H = 3, xi = 5, sigma = 0.002, mobility = 0.03,
        tau_L = 0.12, tau_H = 0.04, gravity = (0.00001, -0.00002)
    )
    phi = [0.5 + 0.3tanh((y - 10 - 0.2cos(2pi * (x - 1) / nx)) / 3) for x in 1:nx, y in 1:ny]
    u = zeros(2, nx, ny)
    u[1, :, :] .= 0.01
    u[2, :, :] .= -0.005
    pressure = fill(0.27, nx, ny)
    for collision in (:bgk, :mrt)
        s = PhaseFieldLBM(phi, p; u, pstar = pressure, collision, threaded = false)
        # Inject an even, zero-mass non-equilibrium stress at every node.
        perturbation = SVector(-0.002, 0.003, -0.002, 0.003, -0.002, 0.001, -0.001, 0.001, -0.001)
        @test sum(perturbation) ≈ 0 atol = 1.0e-17
        for y in 1:ny, x in 1:nx
            s.g[:, x, y] .+= perturbation
        end
        u_old = copy(s.u)
        macroscopic!(s)
        @test s.pstar ≈ pressure atol = 1.0e-14
        @test s.rho ≈ 1 .+ 2s.phi
        @test s.rho != s.pstar
        @test s.tau ≈ p.tau_L .+ s.phi .* (p.tau_H - p.tau_L)
        for (x, y) in ((3, 8), (1, 1), (nx, ny))
            _, geq_old = direct_equilibria(s.lattice, s.phi[x, y], s.pstar[x, y], u_old[:, x, y])
            K = dense_K(s.collision, s.tau[x, y])
            q = K * (s.g[:, x, y] - geq_old) # REGULAR equilibrium
            stress = sum((s.lattice.c[:, i] * s.lattice.c[:, i]') * q[i] for i in 1:9)
            grad = s.grad_phi[:, x, y]
            grad_rho = (p.rho_H - p.rho_L) * grad
            Fmu = -s.tau[x, y] * stress * grad_rho
            chemical = 4 * (12p.sigma / p.xi) * s.phi[x, y] * (s.phi[x, y] - 1) * (s.phi[x, y] - 0.5) -
                (3p.sigma * p.xi / 2) * s.lap_phi[x, y]
            F = chemical * grad + s.rho[x, y] * p.gravity -
                s.pstar[x, y] * Fakhari.CS2 * grad_rho + Fmu
            @test s.force_mu[:, x, y] ≈ Fmu atol = 1.0e-14
            @test s.force[:, x, y] ≈ F atol = 1.0e-14
            @test s.u[:, x, y] ≈ s.lattice.c * s.g[:, x, y] + F / (2s.rho[x, y]) atol = 1.0e-14
        end
        collide!(s)
        x, y = 3, 8
        heq, geq = direct_equilibria(s.lattice, s.phi[x, y], s.pstar[x, y], s.u[:, x, y])
        grad = s.grad_phi[:, x, y]
        normal = grad / norm(grad)
        fp = [(1 - 4(s.phi[x, y] - 0.5)^2) / p.xi * s.lattice.w[i] * dot(s.lattice.c[:, i], normal) for i in 1:9]
        fi = [s.lattice.w[i] * dot(s.lattice.c[:, i], s.force[:, x, y]) / (s.rho[x, y] * Fakhari.CS2) for i in 1:9]
        @test sum(fp) ≈ 0 atol = 1.0e-17
        @test sum(fi) ≈ 0 atol = 1.0e-17
        @test s.h_post[:, x, y] ≈ s.h[:, x, y] -
            (s.h[:, x, y] - (heq - fp / 2)) / (p.mobility / Fakhari.CS2 + 0.5) + fp atol = 1.0e-14
        @test s.g_post[:, x, y] ≈ s.g[:, x, y] -
            dense_K(s.collision, s.tau[x, y]) * (s.g[:, x, y] - (geq - fi / 2)) + fi atol = 1.0e-14
    end
end

@testset "Streaming is a mass-preserving permutation, with RTI walls" begin
    s = PhaseFieldLBM(fill(0.3, 8, 12), Parameters(); threaded = false)
    s.h_post .= reshape(collect(1.0:length(s.h_post)), size(s.h_post))
    s.g_post .= -s.h_post
    mass = sum(s.h_post)
    stream!(s)
    @test sum(s.h) == mass
    @test sort(vec(s.h)) == sort(vec(s.h_post))
    @test s.g == -s.h
    @test s.h[2, 1, 4] == s.h_post[2, 8, 4] # east, periodic seam
    @test s.h[3, 4, 1] == s.h_post[5, 4, 1] # bottom bounce-back
    @test s.h[5, 4, 12] == s.h_post[3, 4, 12] # top bounce-back
    @test s.h[6, 1, 1] == s.h_post[8, 1, 1] # corner diagonal
end

@testset "Uniform states, phase mass and diagnostics" begin
    for T in (Float64, Float32), collision in (:bgk, :mrt), value in (0, 0.37, 1)
        p = Parameters(; T, sigma = 0, gravity = (0, 0))
        phi = fill(T(value), 8, 16)
        s = PhaseFieldLBM(phi, p; collision, pstar = fill(T(0.7), 8, 16), threaded = false)
        step!(s, 100)
        tol = T === Float64 ? 1.0e-12 : 1.0e-5
        @test s.phi ≈ phi atol = tol
        @test maximum(abs, s.u) < tol
        @test s.pstar ≈ fill(T(0.7), 8, 16) atol = tol
        d = diagnostics(s)
        @test d.finite
        @test abs(d.relative_mass_drift) < tol
        @test d.max_mach < tol
        @test d.step == 100
    end
    # Overshoot remains visible, including outside the constructor's input range.
    s = PhaseFieldLBM(fill(1.001, 8, 16), Parameters(); threaded = false)
    @test diagnostics(s).phi_max == 1.001
    s.g[1, 1, 1] = Inf
    @test !diagnostics(s).finite
end

@testset "Viscosity interpolation and RTI dimensionless consistency" begin
    for model in (:linear, :harmonic, :dynamic)
        p = Parameters(tau_L = 0.12, tau_H = 0.04, viscosity_model = model)
        @test relaxation_time(p, 0.0, p.rho_L) ≈ p.tau_L
        @test relaxation_time(p, 1.0, p.rho_H) ≈ p.tau_H
        phi = 0.3
        rho = p.rho_L + phi * (p.rho_H - p.rho_L)
        expected = model === :linear ? 0.12 + phi * (0.04 - 0.12) :
            model === :harmonic ? 1 / (1 / 0.12 + phi * (1 / 0.04 - 1 / 0.12)) :
            ((1 - phi) * p.rho_L * 0.12 + phi * p.rho_H * 0.04) / rho
        @test relaxation_time(p, phi, rho) ≈ expected
    end
    b = rti_parameters(Nx = 256)
    p = b.parameters
    @test (b.Nx, b.Ny) == (256, 1024)
    @test p.xi == 5
    @test p.viscosity_model === :linear
    @test p.rho_H / p.rho_L ≈ 3
    @test (p.rho_H - p.rho_L) / (p.rho_H + p.rho_L) ≈ b.At ≈ 0.5
    @test sqrt(b.L / (b.g * b.At)) ≈ b.t0 ≈ 16000
    @test b.U0 ≈ sqrt(b.g * b.L)
    @test p.rho_H * b.U0 * b.L / b.mu_H ≈ b.Re ≈ 3000
    @test b.mu_H / b.mu_L ≈ b.mu_star ≈ 1
    @test b.mu_H * b.U0 / p.sigma ≈ b.Ca ≈ 0.26
    @test b.U0 * b.L / p.mobility ≈ b.Pe ≈ 1000
    @test p.tau_H * Fakhari.CS2 ≈ b.nu_H
    @test p.tau_L * Fakhari.CS2 ≈ b.nu_L
    @test b.tau_phi * Fakhari.CS2 ≈ p.mobility
    s = initialize_rti(b; threaded = false)
    @test s.collision isa MRTCollision
    @test s.phi[1, 1] < 1.0e-12
    @test s.phi[1, end] > 1 - 1.0e-12
    fronts = interface_fronts(s)
    @test fronts.bubble ≈ 0.1b.L atol = 0.02
    @test fronts.spike ≈ -0.1b.L atol = 0.02
    @test fronts.bubble_crossings == fronts.spike_crossings == 1
    # Actual baseline lattice parameters and half-force initialization.
    @test maximum(abs, s.u) == 0
    step!(s, 10)
    @test diagnostics(s).finite
    @test abs(diagnostics(s).relative_mass_drift) < 1.0e-12
end

@testset "RTI grows under gravity (actual Fig.5 lattice parameters)" begin
    # Moment tests alone cannot detect a solver with no macroscopic evolution.
    # At t*=0.0625 the benchmark has already moved ~0.25 lattice nodes in both
    # directions. Use a threshold comfortably above initialization relaxation.
    setup = rti_parameters(Nx = 256)
    s = initialize_rti(setup)
    initial = interface_fronts(s)
    step!(s, 1000)
    front, d = interface_fronts(s), diagnostics(s)
    @test front.bubble - initial.bubble > 0.1
    @test front.spike - initial.spike < -0.1
    @test d.max_mach > 5.0e-4
    @test d.finite
    @test abs(d.relative_mass_drift) < 1.0e-10
    @test -1.0e-3 < d.phi_min <= d.phi_max < 1.001
end

@testset "Re=6000 tall RTI startup with actual lattice parameters" begin
    setup = rti_parameters(Nx = 256, Re = 6000)
    s = initialize_rti(setup)
    initial = interface_fronts(s)
    # Startup regression only: t*=0.0625, not the 192000-step production run.
    step!(s, 1000)
    front, d = interface_fronts(s), diagnostics(s)
    @test s.collision isa MRTCollision
    @test front.bubble - initial.bubble > 0.1
    @test front.spike - initial.spike < -0.1
    @test d.finite && d.rho_min > 0 && d.tau_min > 0
    @test abs(d.relative_mass_drift) < 1.0e-10
    @test -1.0e-3 < d.phi_min <= d.phi_max < 1.001
    @test 0 < d.max_mach < 0.1
end

@testset "Square RTI preserves vertical resolution and seeds four fingers" begin
    square, baseline = rti_parameters(), rti_parameters(Nx = 256)
    @test (square.Nx, square.Ny) == (1024, 1024)
    @test square.Ny == baseline.Ny
    @test square.L == baseline.L == 256
    @test square.wavelength == 256
    @test square.fingers == 4
    for field in fieldnames(Parameters{Float64})
        @test getfield(square.parameters, field) == getfield(baseline.parameters, field)
    end
    @test square.t0 == baseline.t0
    @test_throws ArgumentError rti_parameters(Nx = 1000)
    @test_throws ArgumentError rti_parameters(Nx = 128)

    # A reduced grid checks all four periods through collision and periodic
    # streaming. A wider domain must evolve like four copies of one period.
    setup = rti_parameters(L = 16)
    repeated = initialize_rti(setup; threaded = false)
    single = initialize_rti(rti_parameters(L = 16, Nx = 16); threaded = false)
    @test size(repeated.phi) == (64, 64)
    @test repeated.phi == repeat(single.phi, 4, 1)
    front = interface_fronts(repeated; wavelength = setup.wavelength)
    @test front.bubble > 0
    @test front.spike < 0 # Nx/2 would wrongly measure another bubble.
    @test_throws ArgumentError interface_fronts(repeated; wavelength = 15)
    step!(repeated, 20)
    step!(single, 20)
    @test repeated.phi == repeat(single.phi, 4, 1)
    @test repeated.g == repeat(single.g, 1, 4, 1)
    @test repeated.u == repeat(single.u, 1, 4, 1)
    @test diagnostics(repeated).finite
    @test abs(diagnostics(repeated).relative_mass_drift) < 1.0e-12
end

@testset "Serial/threaded equivalence" begin
    setup = rti_parameters(L = 16)
    a = initialize_rti(setup; threaded = false)
    b = initialize_rti(setup; threaded = true)
    step!(a, 20)
    step!(b, 20)
    @test a.h == b.h
    @test a.g == b.g
    @test a.phi == b.phi
    @test a.u == b.u
end

module CavityReference
    include("../cavity.jl")
end

@testset "Existing cavity uses the shared lattice unchanged" begin
    lattice = CavityReference.D2Q9Lattice(Float64)
    u = SVector(0.02, -0.01)
    for i in 1:9
        c, w = lattice.c[:, i], lattice.w[i]
        @test CavityReference.f_eq(w, c, u, 1.1) ≈
            w * 1.1 * (1 + 3dot(c, u) + 4.5dot(c, u)^2 - 1.5dot(u, u))
    end
    lbm = CavityReference.D2Q9LBM(
        ones(8, 8), zeros(2, 8, 8),
        CavityReference.D2Q9BGKCollision(1.0)
    )
    CavityReference.macroscopic!(lbm)
    CavityReference.collide!(lbm)
    CavityReference.stream!(lbm)
    @test all(isfinite, lbm.f)
end
