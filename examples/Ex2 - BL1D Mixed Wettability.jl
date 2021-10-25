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

kr_mw = MixedRelPerms(kr_ww, kr_ow)


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

PV = res.A * res.L * res.phi
tmax = PV / bc.value[1] * 3
dt = res.A * grid.dx / bc.value[1] / 50
Int(tmax ÷ dt) + 1

#* Initialize simulation and postprocessing
sim = Simulation(
                1000*ones(grid.Nx)* 6894.76,
                0.2*ones(grid.Nx),
                0.5*ones(grid.Nx),
                tmax,
                dt
)

sim_ww = deepcopy(sim)
sim_ow = deepcopy(sim)

post = FlowResults(sim, grid, tmax, dt)
post_ww = FlowResults(sim, grid, tmax, dt)
post_ow = FlowResults(sim, grid, tmax, dt)


while sim.t <= tmax
    IMPES!(sim, res, kr_mw, fluid, grid, bc, post)
    IMPES!(sim_ww, res, kr_ww, fluid, grid, bc, post_ww)
    IMPES!(sim_ow, res, kr_ow, fluid, grid, bc, post_ow)
end


i = 200
plot(grid.x, post.s[:,i], seriestype=:scatter, label="Mixed-Wet", alpha=0.5)
plot!(grid.x, post_ww.s[:,i], seriestype=:scatter, label="Water-Wet", alpha=0.5)
plot!(grid.x, post_ow.s[:,i], seriestype=:scatter, label="Oil-Wet", alpha=0.5)


#* Recovery Factors
RF = recoveryFactor(grid, post)
RF_ww = recoveryFactor(grid, post_ww)
RF_ow = recoveryFactor(grid, post_ow)
PVs = injPoreVolumes(bc, res, sim)

plot(PVs, RF, label="Mixed Wet", lw=3, alpha=0.5)
plot!(PVs, RF_ww, label="Water Wet", lw=3, alpha=0.5)
plot!(PVs, RF_ow, label="Oil Wet", lw=3, alpha=0.5)