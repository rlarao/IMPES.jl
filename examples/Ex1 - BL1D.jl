using Plots, LinearAlgebra
include("types.jl")
include("flowFunctions.jl")
include("T_matrix.jl")
include("D_matrix.jl")
include("Q_vector.jl")
include("solvers.jl")


#* Define Relative Perms
kr = RelPerms(
        swr=0.2,
        sor=0.2,
        krw0=0.2,
        kro0=1.0,
        nw=3.0,
        no=3.0,
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
    #* Get matrices
    sim.D, sim.d11, sim.d12, sim.d21, sim.d22 = D_matrix(res, fluid, grid, sim)
    sim.T, sim.Tw, sim.To = T_matrix(res, fluid, grid, sim, bc)
    sim.Q, sim.Qw, sim.Qo = Q_vector(res, fluid, grid, wells, sim)
    
    #* Solve
    sim.p = (sim.T + sim.D) \ (sim.D * sim.p + sim.Q)
    sim.s =  sim.s + (sim.d12) \ (- sim.Tw * sim.p + sim.Qw)

    #* Save results
    t += sim.dt
    i += 1

    post.p[:,i] = sim.p
    post.s[:,i] = sim.s
    post.t[i] = t
end

plot(grid.x, post.s[:, end-1], marker=:circle, ylim=(0,1))
plot!(grid.x, post.s[:, 10], marker=:circle, ylim=(0,1))
plot!(grid.x, post.s[:, 1000], marker=:circle, ylim=(0,1))
plot!(grid.x, post.s[:, 500], marker=:circle, ylim=(0,1))

IMPES!(sim, res, fluid, grid, bc, wells)
