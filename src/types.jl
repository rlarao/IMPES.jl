using Base: Float64


mutable struct RelPerms{T <: Float64}
    swr::T
    sor::T
    krw0::T
    kro0::T
    nw::T
    no::T
end

function RelPerms(;swr::T, sor::T, krw0::T, kro0::T, nw::T, no::T) where T <: Float64
    RelPerms{Float64}(swr, sor, krw0, kro0, nw, no)
end


mutable struct Reservoir{T <: Float64}
    L::T
    W::T
    h::T        
    k::T
    phi::T
    kr::RelPerms
    cf::T
end

function Reservoir(;L::T, W::T, h::T, k::T, phi::T, kr::RelPerms, cf::T) where T <: Float64
    Reservoir{Float64}(L, W, h, k, phi, kr, cf)
end


mutable struct Fluids{T <: Real}
    cw::T
    co::T
    μw::T
    μo::T
end
function Fluids(;cw::T, co::T, μw::T, μo::T) where T <: Float64
    Fluids{Float64}(cw, co, μw, μo)
end



mutable struct Grid{T <:Int64, P<:Float64}
    Nx::T
    Ny::T
    dx::P
    dy::P
    x::Vector{P}
    y::Vector{P}
end

function Grid(;res::Reservoir, Nx::T, Ny::T) where T <: Int64
    dx = res.L / Nx 
    dy = res.W / Ny
    x = collect(dx:dx:Nx*dx)
    y = collect(dy:dy:Ny*dy)

    Grid(Nx, Ny, dx, dy, x, y)
end


mutable struct BoundaryConditions
    type:: Vector{Symbol}
    value:: Vector{Real}
end

mutable struct Wells{T<:Vector{Symbol}, P<:Vector{Float64}}
    N::Int64
    type::T
    x::P
    y::P
    rw::P
    constraint::T
    status::T
    label::T
    Q::P
    Pwf::P
    S::P
    i::Vector{Int64}
    j::Vector{Int64}
    l::Vector{Int64}
end

function Wells(;grid::Grid, N::Int64, type::T, x::P, y::P, rw::P,
             constraint::T, label::T, Q::P, Pwf::P) where {T<:Vector{Symbol}, P<:Vector{Float64}}
    Nx = grid.Nx
    dx = grid.dx
    dy = grid.dy

    i = zeros(Int64, N)
    j = zeros(Int64, N)
    l = zeros(Int64, N)

    for k in 1:N
        i[k] = ceil(Int64, x[k] / dx)
        j[k] = ceil(Int64, y[k] / dy)
        l[k] = (j[k] - 1) * Nx + i[k]
    end
    status = [:open]
    S = [0.]
    Wells(N, type, x, y, rw, constraint, status, label, Q, Pwf, S, i, j, l)
end

mutable struct Simulation{T <:Float64}
    p::Vector{T}
    s::Vector{T}
    dt::T
    pc::Vector{T}
    T::Matrix{T}
    Tw::Matrix{T}
    To::Matrix{T}
    Q::Vector{T}
    Qw::Vector{T}
    Qo::Vector{T}
    D::Matrix{T}
    d11::Matrix{T}
    d12::Matrix{T}
    d21::Matrix{T}
    d22::Matrix{T}
end

function Simulation(p::T, s::T, dt::Float64) where T <: Vector{Float64}
    n = length(p)
    
    Simulation(p, s, dt, zeros(n), zeros(n,n), zeros(n,n), zeros(n,n), zeros(n),
                zeros(n), zeros(n), zeros(n,n), zeros(n,n), zeros(n,n), zeros(n,n),
                zeros(n,n))
end


mutable struct PostProcessing
    t::Vector{Int64}
    s::Matrix{Float64}
    p::Matrix{Float64}
   
   function PostProcessing(sim::Simulation, grid::Grid, tmax::Int64, dt::Float64)
    ntimes =  ceil(Int64, tmax / dt)
    p = zeros(grid.Nx*grid.Ny, ntimes)
    s = zeros(grid.Nx*grid.Ny, ntimes)
    
    p[:,1] = sim.p
    s[:,1] = sim.s
    
    t = zeros(ntimes)
    new(t, s, p)
    end

end