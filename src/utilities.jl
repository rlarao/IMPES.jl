function recoveryFactor(grid::Grid, post::FlowResults)
    nsteps = length(post.s[1,:])
    OOIP = sum(1 .- post.s[:,1]) / grid.Nx
    RF = [100 - sum(1 .- post.s[:,i]) / grid.Nx / OOIP * 100  for i in 1:nsteps]
end

function injPoreVolumes(bc::BoundaryConditions, res::Reservoir, sim::Simulation)
    PV = [bc.value[1] * t / (res.L*res.A*res.phi) for t in 0:sim.dt:sim.tmax]
end