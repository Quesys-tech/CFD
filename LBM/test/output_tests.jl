using Test
include("../rayleigh_taylor.jl")

@testset "PVD is available at every completed snapshot" begin
    setup = rti_parameters(L=16)
    s = initialize_rti(setup; threaded=false)
    mktempdir() do dir
        save_snapshot(s, setup, dir)
        path = joinpath(dir,"rti.pvd")
        @test isfile(path)
        @test isfile(joinpath(dir,"step_0.vtr"))
        initial = read(path,String)
        @test occursin("timestep=\"0.0\"",initial)
        @test length(collect(eachmatch(r"<DataSet\s",initial))) == 1

        step!(s,100)
        save_snapshot(s, setup, dir)
        updated = read(path,String)
        @test occursin("step_0.vtr",updated)
        @test occursin("step_100.vtr",updated)
        @test occursin("timestep=\"0.00625\"",updated)
        @test length(collect(eachmatch(r"<DataSet\s",updated))) == 2
        @test read(joinpath(dir,"step_0.vtr")) != read(joinpath(dir,"step_100.vtr"))

        # Restarting in the same directory must not reuse an old time index.
        s0 = initialize_rti(setup; threaded=false)
        save_snapshot(s0, setup, dir)
        restarted = read(path,String)
        @test !occursin("step_100.vtr",restarted)
        @test length(collect(eachmatch(r"<DataSet\s",restarted))) == 1
    end
end

@testset "Runner writes the completed time sequence and motion diagnostics" begin
    mktempdir() do dir
        s = run_rti(L=16,steps=2,output_every=1,diagnostic_every=1,
                    threaded=false,result_dir=dir)
        @test s.iteration == 2
        @test length(collect(eachmatch(r"<DataSet\s",read(joinpath(dir,"rti.pvd"),String)))) == 3
        rows = readlines(joinpath(dir,"diagnostics.csv"))
        @test length(rows) == 4
        @test startswith(rows[end],"2,0.000125,")
    end
end
