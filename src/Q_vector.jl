function Q_vector(res::Reservoir, fluid::Fluids, grid::Grid, well::Wells, sim::Simulation)    
    Qw = zeros(grid.Nx*grid.Ny) 
    Qo = zeros(grid.Nx*grid.Ny) 

    for k in 1:well.N
        l = well.l[k]

        if well.constraint[k] == :Q
            Qw[l] = well.Q[k] # ft3
        end
    end

    Q = Qw + Qo

    return Q, Qw, Qo
end