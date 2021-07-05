

function IMPES!(sim::Simulation, res::Reservoir, fluid::Fluids, grid::Grid,
                bc::BoundaryConditions, wells::Wells)
    #* Get matrices
    sim.D, sim.d11, sim.d12, sim.d21, sim.d22 = D_matrix(res, fluid, grid, sim)
    sim.T, sim.Tw, sim.To = T_matrix(res, fluid, grid, sim, bc)
    sim.Q, sim.Qw, sim.Qo = Q_vector(res, fluid, grid, wells, sim)
    
    #* Solve
    sim.p = (sim.T + sim.D) \ (sim.D * sim.p + sim.Q)
    sim.s =  sim.s + (sim.d12) \ (- sim.Tw * sim.p + sim.Qw)
end