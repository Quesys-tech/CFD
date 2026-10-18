include("fakhari.jl")
using .Fakhari
using WriteVTK
using Printf

"""
Runner presets only: the Fakhari collision, forcing and walls are shared.
high_re_tall matches the square case's lattice updates, not its wall-clock time.
At fixed t0 its final t*=12 includes late-time interaction with the vertical walls.
"""
function rti_case_options(; case = :square, L = 256, t0 = 16000)
    case in (:square, :high_re_tall) ||
        throw(ArgumentError("case must be :square or :high_re_tall"))
    square_steps = round(Int, 3t0)
    if case === :high_re_tall
        return (
            Nx = L, Re = 6000, steps = 4square_steps, output_every = 4000,
            diagnostic_every = 400,
            result_dir = joinpath(@__DIR__, "results", "rti_high_re_tall"),
        )
    end
    return (
        Nx = 4L, Re = 3000, steps = square_steps, output_every = 1000,
        diagnostic_every = 100,
        result_dir = joinpath(@__DIR__, "results", "rti_square"),
    )
end

function write_snapshot(s, setup, result_dir, pvd)
    nx, ny = size(s.phi)
    step = s.iteration
    return vtk_grid(
        joinpath(result_dir, "step_$step"), collect(0.0:(nx - 1)),
        collect((-2setup.L + 0.5):1:(2setup.L - 0.5)), [0.0]
    ) do vtk
        vtk["Phase"] = s.phi
        vtk["Density"] = s.rho
        vtk["NormalizedPressure"] = s.pstar
        vtk["Pressure"] = s.rho .* (Fakhari.CS2 .* s.pstar)
        vtk["Velocity"] = (s.u[1, :, :], s.u[2, :, :], zeros(eltype(s.u), nx, ny))
        vtk["Force"] = (s.force[1, :, :], s.force[2, :, :], zeros(eltype(s.u), nx, ny))
        vtk["ViscousCorrection"] = (
            s.force_mu[1, :, :], s.force_mu[2, :, :],
            zeros(eltype(s.u), nx, ny),
        )
        pvd[step / setup.t0] = vtk
    end
end

function save_snapshot(s, setup, result_dir)
    # vtk_save closes the collection. Reopen its completed index for each
    # snapshot, publishing a reloadable PVD while the simulation is running.
    # Step zero replaces the index from an earlier run in this directory.
    return paraview_collection(joinpath(result_dir, "rti.pvd"); append = s.iteration > 0) do pvd
        write_snapshot(s, setup, result_dir, pvd)
    end
end

function write_diagnostics(io, s, setup)
    d, front = diagnostics(s), interface_fronts(s; wavelength = setup.wavelength)
    println(
        io, join(
            (
                d.step, d.step / setup.t0, d.phase_mass, d.mass_drift,
                d.relative_mass_drift, d.phi_min, d.phi_max, d.max_mach, d.finite,
                d.rho_min, d.tau_min, front.bubble, front.spike,
                front.bubble_crossings, front.spike_crossings,
            ), ','
        )
    )
    flush(io)
    # Report overshoot as measured; never clamp phi or repair distributions.
    if !d.finite || d.rho_min <= 0 || d.tau_min <= 0
        error("Invalid state at step $(d.step): $d")
    end
    return d
end

"""
Run RTI with CPU kernels, horizontal periodicity and vertical no-slip/no-flux
walls. Default: a 1024x1024 square with four wavelengths and Fig.5 material
parameters. Default output includes tstar=0,1,2,3.
"""
function run_rti(;
        case = :square, L = 256, Nx = nothing, Re = nothing, t0 = 16000,
        steps = nothing, output_every = nothing, diagnostic_every = nothing,
        collision = :mrt, threaded = true, T = Float64, result_dir = nothing
    )
    defaults = rti_case_options(; case, L, t0)
    Nx = something(Nx, defaults.Nx)
    Re = something(Re, defaults.Re)
    steps = something(steps, defaults.steps)
    output_every = something(output_every, defaults.output_every)
    diagnostic_every = something(diagnostic_every, defaults.diagnostic_every)
    result_dir = something(result_dir, defaults.result_dir)
    steps >= 0 && output_every > 0 && diagnostic_every > 0 ||
        throw(ArgumentError("invalid steps or output intervals"))
    setup = rti_parameters(; T, L, Nx, Re, t0)
    s = initialize_rti(setup; collision, threaded)
    initial_front = interface_fronts(s; wavelength = setup.wavelength)
    mkpath(result_dir)
    @printf(
        "RTI: %d x %d, %d steps, final t*=%.4f (t0=%.0f steps)\n",
        setup.Nx, setup.Ny, steps, steps / setup.t0, setup.t0
    )
    @printf(
        "reference L=%.0f, wavelength=%d, initial fingers=%d\n",
        setup.L, setup.wavelength, setup.fingers
    )
    @printf("case=%s Re=%.0f lattice updates=%d\n", string(case), setup.Re, setup.Nx * setup.Ny * steps)
    if steps / setup.t0 > 3
        println("Extended run beyond t*=3: vertical walls can affect the late-time evolution.")
    end
    if steps / setup.t0 < 0.1
        println(
            "Short startup run: the interface displacement can be subpixel. " *
                "Use --steps=$(round(Int, 3setup.t0)) to inspect t*=0,1,2,3."
        )
    end
    # Persist the exact derived lattice parameters alongside every run.
    open(joinpath(result_dir, "parameters.txt"), "w") do io
        println(io, "Fakhari et al. (2017), velocity-based D2Q9 phase-field LBM")
        println(io, "dx=dt=c=1, cs2=1/3; phi_L=0, phi_H=1")
        println(io, "case = $case")
        println(io, "steps = $steps")
        println(io, "final_tstar = $(steps / setup.t0)")
        println(io, "lattice_updates = $(setup.Nx * setup.Ny * steps)")
        println(io, "output_every = $output_every, diagnostic_every = $diagnostic_every")
        println(io, "collision=$collision, viscosity_model=$(s.parameters.viscosity_model)")
        println(io, "initial signed-distance approximation: d=y-yI(x)")
        println(io, "pressure: columnwise hydrostatic, zero at top wall")
        println(io, "stress equilibrium uses previous-step velocity (explicit)")
        for (key, value) in pairs(setup)
            key === :parameters || println(io, "$key = $value")
        end
        println(io, "tau_L = $(s.parameters.tau_L)")
        println(io, "tau_H = $(s.parameters.tau_H)")
        println(io, "xi = $(s.parameters.xi)")
    end
    snapshots = Set(round.(Int, (0, setup.t0, 2setup.t0, 3setup.t0)))
    open(joinpath(result_dir, "diagnostics.csv"), "w") do io
        println(
            io, "step,tstar,phase_mass,mass_drift,relative_mass_drift," *
                "phi_min,phi_max,max_mach,finite,rho_min,tau_min,bubble_y,spike_y," *
                "bubble_crossings,spike_crossings"
        )
        for iteration in 0:steps
            output = iteration % output_every == 0 || iteration in snapshots || iteration == steps
            if iteration % diagnostic_every == 0 || output
                d = write_diagnostics(io, s, setup)
                if output
                    front = interface_fronts(s; wavelength = setup.wavelength)
                    @printf(
                        "step=%d t*=%.4f mass drift=%+.3e phi=[%.6f,%.6f] Ma=%.4g bubble=%.3f (delta=%+.3f) spike=%.3f (delta=%+.3f)\n",
                        iteration, iteration / setup.t0, d.relative_mass_drift,
                        d.phi_min, d.phi_max, d.max_mach, front.bubble,
                        front.bubble - initial_front.bubble, front.spike,
                        front.spike - initial_front.spike
                    )
                end
            end
            output && save_snapshot(s, setup, result_dir)
            iteration == steps || step!(s)
        end
    end
    return s
end

function main(args = ARGS)
    if "--help" in args
        println(
            """
            julia --project=. --threads=auto LBM/rayleigh_taylor.jl [options]
              --case=square           square or high_re_tall
              --nx=N                  horizontal nodes; square: 4L; high_re_tall: L
              --length=256            reference L and wavelength; Ny=4L (even integer)
              --re=3000               Reynolds number; high_re_tall default: 6000
              --t0=16000              reference time in lattice steps
              --steps=48000           square: round(3*t0); high_re_tall: 4*round(3*t0)
              --output-every=1000     high_re_tall: 4000 (+ t*=0,1,2,3 and final)
              --diagnostic-every=100  high_re_tall: 400; mass, extrema, Mach, fronts
              --collision=mrt        mrt or bgk
              --output-dir=PATH       case-specific directory under LBM/results
              --serial               disable threading
            """
        )
        return
    end
    options = Dict{Symbol, Any}()
    for arg in args
        if arg == "--serial"
            options[:threaded] = false
            continue
        end
        parts = split(arg, '='; limit = 2)
        length(parts) == 2 || throw(ArgumentError("expected --option=value: $arg"))
        key, value = parts
        if key == "--case"
            options[:case] = Symbol(value)
        elseif key == "--re"
            options[:Re] = parse(Float64, value)
        elseif key == "--nx"
            options[:Nx] = parse(Int, value)
        elseif key == "--length"
            options[:L] = parse(Int, value)
        elseif key == "--t0"
            options[:t0] = parse(Float64, value)
        elseif key == "--steps"
            options[:steps] = parse(Int, value)
        elseif key == "--output-every"
            options[:output_every] = parse(Int, value)
        elseif key == "--diagnostic-every"
            options[:diagnostic_every] = parse(Int, value)
        elseif key == "--collision"
            options[:collision] = Symbol(value)
        elseif key == "--output-dir"
            options[:result_dir] = abspath(value)
        else
            throw(ArgumentError("unknown option: $key"))
        end
    end
    return run_rti(; options...)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
