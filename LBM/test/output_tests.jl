using Test
include("../rayleigh_taylor.jl")

@testset "High-Re tall preset has the square case's lattice-update budget" begin
    square = rti_case_options()
    tall = rti_case_options(case = :high_re_tall)
    a = rti_parameters(Nx = square.Nx, Re = square.Re)
    b = rti_parameters(Nx = tall.Nx, Re = tall.Re)
    @test (b.Nx, b.Ny) == (256, 1024)
    @test b.Re == 6000
    @test tall.steps == 192000
    @test tall.steps / b.t0 == 12
    @test a.Nx * a.Ny * square.steps == b.Nx * b.Ny * tall.steps == 50331648000
    @test tall.steps / tall.output_every == square.steps / square.output_every
    @test tall.steps / tall.diagnostic_every == square.steps / square.diagnostic_every
    @test b.U0 == a.U0 && b.g == a.g && b.t0 == a.t0
    @test b.parameters.viscosity_model === :linear
    @test b.parameters.tau_H ≈ a.parameters.tau_H / 2
    @test b.parameters.tau_L ≈ a.parameters.tau_L / 2
    @test b.sigma ≈ a.sigma / 2
    @test b.mobility == a.mobility
    @test b.parameters.rho_H * b.U0 * b.L / b.mu_H ≈ 6000
    @test b.mu_H * b.U0 / b.sigma ≈ 0.26
    @test b.U0 * b.L / b.mobility ≈ 1000
    @test_throws ArgumentError rti_case_options(case = :unknown)
    @test_throws ArgumentError rti_parameters(Re = 0)
    @test_throws ArgumentError rti_parameters(Re = Inf)
    # Case defaults scale the work budget with L and t0 without changing Ny=4L.
    tiny = rti_case_options(case = :high_re_tall, L = 16, t0 = 100)
    @test tiny.Nx == 16 && tiny.steps == 1200
end

@testset "High-Re CLI overrides and VTK geometry (two small-grid steps)" begin
    mktempdir() do dir
        s = main(
            [
                "--re=7500", "--length=16", "--steps=2", "--output-every=1",
                "--diagnostic-every=1", "--case=high_re_tall", "--serial",
                "--output-dir=$dir",
            ]
        )
        @test size(s.phi) == (16, 64)
        @test s.collision isa MRTCollision
        @test s.iteration == 2 && diagnostics(s).finite
        params = read(joinpath(dir, "parameters.txt"), String)
        @test occursin("Re = 7500.0", params)
        @test occursin("case = high_re_tall", params)
        @test occursin("lattice_updates = 2048", params)
        xml = first(split(String(read(joinpath(dir, "step_2.vtr"))), "<AppendedData"; limit = 2))
        @test occursin("WholeExtent=\"0 15 0 63 0 0\"", xml)
        @test length(collect(eachmatch(r"<DataSet\s", read(joinpath(dir, "rti.pvd"), String)))) == 3
    end
end

@testset "PVD is available at every completed snapshot" begin
    setup = rti_parameters(L = 16)
    s = initialize_rti(setup; threaded = false)
    mktempdir() do dir
        save_snapshot(s, setup, dir)
        path = joinpath(dir, "rti.pvd")
        @test isfile(path)
        @test isfile(joinpath(dir, "step_0.vtr"))
        initial = read(path, String)
        @test occursin("timestep=\"0.0\"", initial)
        @test length(collect(eachmatch(r"<DataSet\s", initial))) == 1

        step!(s, 100)
        save_snapshot(s, setup, dir)
        updated = read(path, String)
        @test occursin("step_0.vtr", updated)
        @test occursin("step_100.vtr", updated)
        @test occursin("timestep=\"0.00625\"", updated)
        @test length(collect(eachmatch(r"<DataSet\s", updated))) == 2
        @test read(joinpath(dir, "step_0.vtr")) != read(joinpath(dir, "step_100.vtr"))

        # Restarting in the same directory must not reuse an old time index.
        s0 = initialize_rti(setup; threaded = false)
        save_snapshot(s0, setup, dir)
        restarted = read(path, String)
        @test !occursin("step_100.vtr", restarted)
        @test length(collect(eachmatch(r"<DataSet\s", restarted))) == 1
    end
end

@testset "Runner writes the completed time sequence and motion diagnostics" begin
    mktempdir() do dir
        s = run_rti(
            L = 16, steps = 2, output_every = 1, diagnostic_every = 1,
            threaded = false, result_dir = dir
        )
        @test s.iteration == 2
        @test size(s.phi) == (64, 64)
        @test length(collect(eachmatch(r"<DataSet\s", read(joinpath(dir, "rti.pvd"), String)))) == 3
        rows = readlines(joinpath(dir, "diagnostics.csv"))
        @test length(rows) == 4
        @test startswith(rows[end], "2,0.000125,")
        values = split(rows[2], ',')
        @test parse(Float64, values[12]) > 0
        @test parse(Float64, values[13]) < 0
        parameters = read(joinpath(dir, "parameters.txt"), String)
        @test occursin("Nx = 64", parameters)
        @test occursin("Ny = 64", parameters)
        @test occursin("fingers = 4", parameters)
        # Coordinates and point counts must describe a square, not a stretched
        # display of the old tall mesh. Inspect only the XML before raw data.
        vtk = String(read(joinpath(dir, "step_0.vtr")))
        xml = first(split(vtk, "<AppendedData"; limit = 2))
        @test occursin("WholeExtent=\"0 63 0 63 0 0\"", xml)
    end
end
