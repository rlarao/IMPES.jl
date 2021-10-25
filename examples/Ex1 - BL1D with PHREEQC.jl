# using Revise
using IMPES
using JPhreeqc
# using Plots, LinearAlgebra, IMPES
# using AnalyticalEOR

#* Define Relative Perms
kr = RelPerms(
        swr= 0.2,
        sor= 0.2,
        krw0= 0.3,
        kro0= 1.0,
        nw= 2.0,
        no= 3.0,
            )

#* Define reservoir
res = Reservoir(
    L = 300.,
    W = 200.0 * 0.3048,
    h = 50.0 * 0.3048,
    k = 100.0 * 9.869233e-16,
    phi = 0.2,
    kr = kr,
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
                tmax,
                dt
)

post = FlowResults(sim, grid, tmax, dt)


ncomps = 2
cinj = ones(ncomps) * 1.
ci = zeros(ncomps, grid.Nx)
ntimes = ceil(Int64, sim.tmax / sim.dt)
c = zeros(ncomps, grid.Nx, ntimes+1)
c[:,:, 1] = ci
ctrans = CompTransport(c, ncomps, cinj)




while sim.t < tmax
    IMPES!(sim, res, fluid, grid, bc, post)
    IMPEC!(ctrans, sim, res, grid, post)
end


using Plots

i = 1000
plot(grid.x, post.s[:,i], seriestype=:scatter)
plot(grid.x, ctrans.c[1,:,i], seriestype=:scatter)
plot!(grid.x, ctrans.c[2,:,i], seriestype=:scatter)



