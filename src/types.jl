using Base: Float64

struct RelPerms{T <: Float64}
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


struct Reservoir{T <: Float64}
    L::T
    W::T
    h::T        
    k::T
    phi::T
    kr::RelPerms
end

function Reservoir(;L::T, W::T, h::T, k::T, phi::T, kr::RelPerms) where T <: Float64
    Reservoir{Float64}(L, W, h, k, phi, kr)
end


struct Fluids{T <: Real}
    μw::T
    μo::T
end
function Fluids(;μw::T, μo::T) where T <: Float64
    Fluids{Float64}(μw, μo)
end



mutable struct Grid{T <:Int64, P<:Float64}
    Nx::T
    dx::P
    x::Vector{P}
end

function Grid(;res::Reservoir, Nx::T) where T <: Int64
    dx = res.L / Nx 
    x = collect(dx/2:dx:res.L-dx/2)

    Grid(Nx, dx, x)
end


struct BoundaryConditions
    type:: Vector{Symbol}
    value:: Vector{Real}
end


mutable struct Simulation{T <:Float64}
    p::Vector{T}
    s::Vector{T}
    v::Vector{T}
    i::Int64
    t::T
    tmax::T
    dt::T
end

function Simulation(p::T, s::T, tmax, dt::Float64) where T <: Vector{Float64}
    n = length(p)
    
    Simulation(p, s, zeros(n+1), 1, 0., tmax, dt)
end


mutable struct FlowResults
    t::Vector{Float64}
    s::Matrix{Float64}
    p::Matrix{Float64}

   function FlowResults(sim::Simulation, grid::Grid, tmax::Float64, dt::Float64)
    ntimes =  ceil(Int64, tmax / dt)
    p = zeros(grid.Nx, ntimes+1)
    s = zeros(grid.Nx, ntimes+1)
    
    p[:,1] = sim.p
    s[:,1] = sim.s
    t = collect(0:dt:ntimes*dt)
    
    new(t, s, p)
    end

end


# mutable struct CompTransport
#     c::Matrix{Float64}
#     ncomps::Int64
#     cinj::Vector{Float64}

#     function CompTransport(ncomps::Int64, cinj::Vector{Float64}, ci::Matrix{Float64}, sim::Simulation)
#         n = length(sim.p)
#         ntimes = ceil(Int64, sim.tmax / sim.dt)
#         c = zeros(ncomps, n, ntimes+1)
#         c[:,:, 1] = ci

#         new(c, ncomps, cinj)
#     end
# end


mutable struct CompTransport{N}
    c:: Array{Float64, N}
    ncomps::Int64
    cinj::Vector{Float64}

end
