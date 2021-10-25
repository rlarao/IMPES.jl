using Base: Float64

#* Flow functions
abstract type RelPerms end

mutable struct SingleRelPerms <: RelPerms
    swr::Float64
    sor::Float64
    krw0::Float64
    kro0::Float64
    nw::Float64
    no::Float64
end

mutable struct MixedRelPerms <: RelPerms
    ww::SingleRelPerms
    ow::SingleRelPerms
end

function SingleRelPerms(;swr::T, sor::T, krw0::T, kro0::T, nw::T, no::T) where T <: Float64
    SingleRelPerms(swr, sor, krw0, kro0, nw, no)
end


#* Reservoir
abstract type Reservoir end

struct Cylinder <: Reservoir
    L   ::Float64
    D   ::Float64
    A   ::Float64        
    V   ::Float64        
    k   ::Float64
    phi ::Float64
end

struct Block <: Reservoir
    L   ::Float64
    W   ::Float64
    h   ::Float64
    A   ::Float64        
    V   ::Float64
    k   ::Float64
    phi ::Float64
end

function Reservoir(L::Float64, D::Float64, k::Float64, phi::Float64)
    A = π * D ^2 /4
    V = A * L
    Cylinder(L, D, A, V, k, phi)
end

function Reservoir(L::Float64, W::Float64, h::Float64, k::Float64, phi::Float64)
    V = L * W * h
    A = W * h
    Block(L, W, h, A, V, k, phi)
end

#* Fluid
struct Fluids{T <: Real}
    μw::T
    μo::T
end
function Fluids(;μw::T, μo::T) where T <: Float64
    Fluids{Float64}(μw, μo)
end

#* Grid
abstract type Grid end

mutable struct Grid1D <: Grid
    Nx  ::Int64
    N   ::Int64
    dx  ::Float64
    x   ::Vector{Float64}
end

mutable struct Grid2D <: Grid
    Nx  ::Int64
    Ny  ::Int64
    N   ::Int64
    dx  ::Float64
    dy  ::Float64
    x   ::Vector{Float64}
    y   ::Vector{Float64}
end

function Grid(res::Reservoir, Nx::Int64)
    dx = res.L / Nx
    x = collect(dx/2:dx:res.L-dx/2)

    Grid1D(Nx, Nx, dx, x)
end

function Grid(res::Reservoir, Nx::Int64, Ny::Int64)
    dx = res.L / Nx
    dy = res.W / Ny
    x = collect(dx/2:dx:res.L-dx/2)
    y = collect(dy/2:dy:res.L-dy/2)
    N = Nx * Ny

    Grid2D(Nx, Ny, N, dx, dy, x, y)
end

# function Grid(;res::Reservoir, Nx::T, Ny) where T <: Int64
#     dx = res.L / Nx 
#     x = collect(dx/2:dx:res.L-dx/2)

#     Grid(Nx, dx, x)
# end

struct BoundaryConditions
    type:: Vector{Symbol}
    value:: Vector{Real}
end

abstract type Simulation end

mutable struct Simulation1D <: Simulation
    p   ::Vector{Float64}
    s   ::Vector{Float64}
    θ   ::Vector{Float64}
    v   ::Vector{Float64}
    i   ::Int64
    t   ::Float64
    tmax::Float64
    dt  ::Float64
end

mutable struct Simulation2D <: Simulation
    p   ::Matrix{Float64}
    s   ::Matrix{Float64}
    θ   ::Matrix{Float64}
    v   ::Matrix{Float64}
    i   ::Int64
    t   ::Float64
    tmax::Float64
    dt  ::Float64
end

function Simulation(p::T, s::T, θ::T, tmax::Float64, dt::Float64) where T <: Vector{Float64}
    n = length(p)
    Simulation1D(p, s, θ, zeros(n+1), 1, 0., tmax, dt)
end


function Simulation(p::T, s::T, θ::T, tmax::Float64, dt::Float64) where T <: Matrix{Float64}
    Nx, Ny = size(p)
    Simulation2D(p, s, θ, zeros(Nx, Ny), 1, 0., tmax, dt)
end

mutable struct FlowResults
    t::Vector{Float64}
    p::Matrix{Float64}
    s::Matrix{Float64}
    θ::Matrix{Float64} 

   function FlowResults(sim::Simulation, grid::Grid, tmax::Float64, dt::Float64)
    ntimes =  Int(tmax ÷ dt)
    p = zeros(grid.Nx, ntimes+1)
    s = zeros(grid.Nx, ntimes+1)
    θ = zeros(grid.Nx, ntimes+1)

    t = collect(0:dt:ntimes*dt)
    p[:,1] = sim.p
    s[:,1] = sim.s
    θ[:,1] = sim.θ
    
    new(t, p, s, θ)
    end
end

mutable struct CompTransport{N}
    c:: Array{Float64, N}
    ncomps::Int64
    cinj::Vector{Float64}
end
