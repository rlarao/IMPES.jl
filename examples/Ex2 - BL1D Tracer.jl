using Plots, LinearAlgebra, IMPES
using AnalyticalEOR


swr=0.2
sor=0.2
krw0=0.3
kro0=1.0
nw=2.0
no =3.0

#* Define Relative Perms
kr = IMPES.RelPerms(
        swr= swr,
        sor= sor,
        krw0= krw0,
        kro0= kro0,
        nw= nw,
        no= no,
            )

#* Define reservoir
res = Reservoir(
    L = 1000.0,
    W = 200.0,
    h = 50.0,
    k = 100.0,
    phi = 0.2,
    kr = kr,
    cf = 0.0,
    )

#* Define fluids
fluid = Fluids(
            cw=1e-5,
            co=1e-5,
            μw=1.0,
            μo=1.0,
)

#* Define grid
grid = Grid(
            res=res,
            Nx = 100,
            Ny = 1
            )

#* Define BC
bc = BoundaryConditions(
                        [:neumann, :dirichlet, :neumann, :neumann],
                        [0.0, 100.0, 0.0, 0.0]
                        )

#* Define wells
wells = Wells(
            grid=grid,
            N = 1,
            type = [:V],
            x = [5.0],
            y = [100.0],
            rw = [0.5],
            constraint = [:Q],
            label = [:inj],
            Q = [426.5],
            Pwf = [100.0],
)

tmax = 4250 # days
dt = 1.0

#* Initialize simulation and posprocessing
sim = Simulation(
                1000*ones(grid.Nx * grid.Ny),
                0.2*ones(grid.Nx * grid.Ny), 
                dt
)

post = PostProcessing(
    sim,
    grid, 
    tmax,
    dt
)


t = 0
i = 0
while t < tmax
    IMPES!(sim, res, fluid, grid, bc, wells)

    #* Save results
    t += sim.dt
    i += 1

    post.p[:,i] = sim.p
    post.s[:,i] = sim.s
    post.t[i] = t
end



#* Analytical solution
si = 0.2
sj = 0.8
μo = 1e-3
μw = 1e-3

kr2 = AnalyticalEOR.RelPerms(
                            swr= swr,
                            sor= sor,
                            krw0= krw0,
                            kro0= kro0,
                            nw= nw,
                            no= no,
                                )


wf = solve_waterflooding(si, sj, kr2, μw, μo)


plot_fw(wf)
V = res.h * res.W * res.L
Q = 426.5

td = post.t * Q / V / res.phi 
td
i = 1500



anim2 = @animate for i = 1:10:4250
    plot_sw_profile(wf, td[i])
    plot!(grid.x / res.L, post.s[:, i], seriestype=:scatter, alpha=0.6, label="IMPES")
end

gif(anim, "1DBL_comparison.gif")
