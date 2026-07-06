using StaticArrays
using WriteVTK
using LinearAlgebra
abstract type AbstractLattice{D,Q,T} end
abstract type AbstractCollision{D,Q,T} end

struct D2Q9Lattice{T} <: AbstractLattice{2,9,T}
    w::SVector{9,T}
    c::SMatrix{2,9,T,18}
    function D2Q9Lattice(::Type{T}) where {T}
        w = @SVector [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36]
        c = @SMatrix [0 1 0 -1 0 1 -1 -1 1;
            0 0 1 0 -1 1 1 -1 -1]
        new{T}(w, c)
    end
end

"""
    f_eq(w_k, c_k, u, rho)

Compute the equilibrium distribution function.

# Arguments
- `w_k`: Weight for the k-th direction.
- `c_k`: Discrete velocity vector for the k-th direction.
- `u`: Macroscopic velocity vector.
- `rho`: Macroscopic density.
"""
function f_eq(w_k, c_k, u, rho)
    w_k * rho * (1 + 3c_k ⋅ u + 9 / 2 * (c_k ⋅ u)^2 - 3 / 2 * (u ⋅ u))
end

struct D2Q9BGKCollision{T} <: AbstractCollision{2,9,T}
    omega::T
    function D2Q9BGKCollision(omega::T) where {T}
        new{T}(omega)
    end
end

struct D2Q9LBM{T,A2D<:AbstractArray{T,2},
    A3D<:AbstractArray{T,3}}
    lattice::D2Q9Lattice{T}
    collision::D2Q9BGKCollision{T}
    f::A3D # M*N*9 
    f_post::A3D# M*N*9 
    rho_u::A3D# M*N*2
    rho::A2D# M*N
    u::A3D # M*N*2
    function D2Q9LBM(rho::A2D, u::A3D, collision::D2Q9BGKCollision{T}) where {T,A2D<:AbstractArray{T,2},A3D<:AbstractArray{T,3}}
        N_x, N_y = size(rho)
        @assert N_x == size(u, 2) && N_y == size(u, 3) "rho and u must have the same dimensions"
        lattice = D2Q9Lattice(T)
        f = Array{T}(undef, 9, N_x, N_y)
        for jy in 1:N_y, jx in 1:N_x, i in 1:9
            f[i, jx, jy] = f_eq(lattice.w[i], lattice.c[:, i], SA{T}[u[1, jx, jy], u[2, jx, jy]], rho[jx, jy])
        end
        f_post = zeros(T, 9, N_x, N_y)
        rho_u = zeros(T, 9, N_x, N_y)
        new{T,A2D,A3D}(lattice, collision, f, f_post, rho_u, rho, u)
    end
end

function macroscopic!(lbm::D2Q9LBM{T}) where {T}
    for jy in 1:size(lbm.rho, 2), jx in 1:size(lbm.rho, 1)
        lbm.rho[jx, jy] = 0.0
        lbm.u[1, jx, jy] = 0.0
        lbm.u[2, jx, jy] = 0.0
        for i in 1:9
            c_k = lbm.lattice.c[:, i]
            lbm.rho[jx, jy] += lbm.f[i, jx, jy]
            @. lbm.u[:, jx, jy] += lbm.f[i, jx, jy] * c_k
        end
        lbm.u[:, jx, jy] /= lbm.rho[jx, jy]
    end
end

function collide!(lbm::D2Q9LBM{T}) where {T}
    for jy in 1:size(lbm.rho, 2), jx in 1:size(lbm.rho, 1)
        rho = lbm.rho[jx, jy]
        u = @SVector [lbm.u[1, jx, jy], lbm.u[2, jx, jy]]
        for i in 1:9
            f_eq_k = f_eq(lbm.lattice.w[i], lbm.lattice.c[:, i], u, rho)
            lbm.f_post[i, jx, jy] = lbm.f[i, jx, jy] - lbm.collision.omega * (lbm.f[i, jx, jy] - f_eq_k)
        end
    end
end

function stream!(lbm::D2Q9LBM{T}) where {T}
    N_x, N_y = size(lbm.rho)
    k_opp = [1, 4, 5, 2, 3, 8, 9, 6, 7]
    for jy in 1:N_y, jx in 1:N_x, i in 1:9
        c_i = lbm.lattice.c[:, i]
        jx_pre = jx - Int(c_i[1])
        jy_pre = jy - Int(c_i[2])

        if jx_pre < 1 || jx_pre > N_x || jy_pre < 1
            lbm.f[i, jx, jy] = lbm.f_post[k_opp[i], jx, jy]
        elseif jy_pre > N_y
            U_in = SA{T}[0.05, 0.0]
            rho_in = one(T)  # or lbm.rho[jx, jy], depending on your model
            lbm.f[i, jx, jy] = lbm.f_post[k_opp[i], jx, jy] + 6 * lbm.lattice.w[i] * rho_in * (c_i ⋅ U_in)
        else
            lbm.f[i, jx, jy] = lbm.f_post[i, jx_pre, jy_pre]
        end
    end
end

function main()
    U_max = 0.05
    Re = 100.0
    Nx = 100
    Ny = Nx
    nu = U_max*Nx / Re
    tau = Float32(3nu + 0.5)
    @show tau
    omega = 1.0 / tau
    collision = D2Q9BGKCollision(Float32(omega))
    rho = ones(Float32, Nx, Ny)
    u = zeros(Float32, 2, Nx, Ny)
    lbm = D2Q9LBM(rho, u, collision)
    result_dir = joinpath(@__DIR__, "results/2dtaylor_green")
    mkpath(result_dir)
    paraview_collection(joinpath(result_dir, "Nx$(Nx)_Ny$(Ny).pvd")) do pvd
        for t in 0:10000
            macroscopic!(lbm)
            if t % 100 == 0
                vtk_dst = joinpath(result_dir, "step_$t.vtr")
                @show vtk_dst
                vtk_grid(vtk_dst, collect(0.5:1:(Nx-0.5)), collect(0.5:1:(Ny-0.5)), [0.0]) do vtk
                    vtk["Density"] = lbm.rho
                    vtk["Velocity"] = (lbm.u[1, :, :], lbm.u[2, :, :], zeros(eltype(lbm.u), Nx, Ny))
                    vtk["Mass"] = Nx*Ny*sum(lbm.rho)
                    vtk["KineticEnergy"] = 0.5*sum(lbm.rho .* (lbm.u[1, :, :] .^ 2 + lbm.u[2, :, :] .^ 2))
                    pvd[t] = vtk
                end
            end
            collide!(lbm)
            stream!(lbm)
        end
    end
end

main()
