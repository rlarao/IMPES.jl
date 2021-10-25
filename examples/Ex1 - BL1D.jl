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
    L = 300.,
    A = 200 * 50 * 0.3048^2,
    k = 100.0 * 9.869233e-16,
    phi = 0.2,
    )

#* Define fluids
fluid = Fluids(
            μw=1e-3,
            μo=5e-3,
)

#* Define grid
grid = Grid(
            res=res,
            Nx = 100,
            )

#* Define BC
bc = BoundaryConditions(
                        [:neumann, :dirichlet],
                        [426.5 * 0.3048 ^ 3 / 3600 / 24, 1000 * 6894.76]
                        )

tmax = 4000. * 3600 * 24 # days
dt = 1. * 3600 * 24

#* Initialize simulation and postprocessing
sim = Simulation(
                1000*ones(grid.Nx)* 6894.76,
                0.2*ones(grid.Nx),
                1*ones(grid.Nx),
                tmax,
                dt
)

post = FlowResults(sim, grid, tmax, dt)


while sim.t < tmax
    IMPES!(sim, res, kr_ww, fluid, grid, bc, post)
end


i = 1000
plot(grid.x, post.s[:,i], seriestype=:scatter)
plot(grid.x, post.p[:,i], seriestype=:scatter)



