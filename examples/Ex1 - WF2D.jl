using Revise
using IMPES
using Plots

# using Plots, LinearAlgebra, IMPES
# using AnalyticalEOR

#* Define Relative Perms
kr_ww = SingleRelPerms(
        swr= 0.11,
        sor= 0.10,
        krw0= 0.2,
        kro0= 0.8,
        nw= 3.0,
        no= 2.0,
            )

kr_ow = SingleRelPerms(
        swr= 0.11,
        sor= 0.16,
        krw0= 0.4,
        kro0= 0.5,
        nw= 2.0,
        no= 3.0,
        )

#* Define reservoir
res = Reservoir(
    300.,
    200 * 0.3048,
    50 * 0.3048,
    100.0 * 9.869233e-16,
    0.2,
    )

#* Define fluids
fluid = Fluids(
            μw=1e-3,
            μo=5e-3,
)

#* Define grid
grid = Grid(
            res,
            100,
            )


grid = Grid(res,
                10,
                10)
#* Define BC


bc = BoundaryConditions(
                        [:neumann, :neumann, :neumann, :neumann],
                        [0., 0., 0., 0.]
                        )

sim = Simulation(zeros(grid.Nx, grid.Ny),
                0.2 * ones(grid.Nx, grid.Ny),
                zeros(grid.Nx, grid.Ny),
                100.,
                1.)



function Thalf(l::Int64, m::Int64, res::Reservoir,
            kr::RelPerms, fluid::Fluids, grid::Grid, sim::Simulation)
    Nx = grid.Nx
    Δx = grid.dx
    Δy = grid.dy
    k = res.k
    h = res.h
    μw = fluid.μw
    μo = fluid.μo
    s = sim.s
    p = sim.p         
    θ = sim.θ

    # lj = ceil(Int64,l / Nx)
    # li = l - (lj-1) * Nx

    # mj = ceil(Int64, m / Nx)
    # mi = m - (mj-1) * Nx

    lj= l ÷ Nx + 1
    li= l - (lj - 1) * Nx
    
    mj= l ÷ Nx + 1
    mi= m - (mj - 1) * Nx

    if p[l] >= p[m]
        krw, kro = kr(s[l], θ)
    else
        krw, kro = kr(s[m], θ)
    end

    if lj==mj # Get x-direction transmissibilities
        Tw, To = (krw/μw, kro/μo ) .* (k * Δy * h / Δx)
    elseif li==mi # Get y-direction transmissibilities
        Tw, To = (krw/μw, kro/μo ) .* (k * Δx * h / Δy)
    elseif lj==mj && li==mi
        if l-1 % Nx == 0 || l % Nx == 0
            Tw, To = (krw/μw, kro/μo) .* (k * Δy * h / Δx)
        else
            Tw, To = (krw/μw, kro/μo) .* (k * Δx * h / Δy)
        end
    end
end

T(l,m) = Thalf(l,m, res, kr_ww, fluid, grid, sim)

Nx = grid.Nx
Ny = grid.Ny
N = grid.N

Tw = zeros(N, N)
To = zeros(N, N)

for l in 1:N

# l = 11

    if (l-1) % Nx != 0 #* Not in Left
        Tw[l, l-1], To[l, l-1] = .-T(l,l-1)
        Tw[l,l] -= Tw[l,l-1]
        To[l,l] -= To[l,l-1]
    end

    if l % Nx != 0 #* Not in Right
        Tw[l, l+1], To[l, l+1] = .-T(l,l+1)
        Tw[l,l] -= Tw[l,l+1]
        To[l,l] -= To[l,l+1]
    end

    if ceil(Int64, l / Nx) != 1 #* Not in Bottom
        Tw[l, l-Nx], To[l, l-Nx] = .-T(l,l-Nx)
        Tw[l,l] -= Tw[l,l-Nx]
        To[l,l] -= To[l,l-Nx]
    end

    if ceil(Int64, l / Nx) != Ny #* Not in Top
        Tw[l, l+Nx], To[l, l+Nx] = .-T(l,l+Nx)
        Tw[l,l] -= Tw[l,l+Nx]
        To[l,l] -= To[l,l+Nx]
    end
end


Tt = Tw .+ To
sparse(Tt)

Q = zeros(grid.N)
Q[1] = 426.5 * 0.3048 ^ 2 /3600 / 24
Q[end] = -426.5 * 0.3048 ^ 2 /3600 / 24

Q *= grid.dx * grid.dy

p = Tt \ Q
plot(reshape(p, grid.Nx, grid.Ny) / 6894, linetype=:contourf)

p
