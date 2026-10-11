using Test
include("../rayleigh_taylor.jl")

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
