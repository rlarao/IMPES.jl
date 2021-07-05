function D_matrix(res::Reservoir, fluid::Fluids, grid::Grid, sim::Simulation)
    Nx = grid.Nx
    Ny = grid.Ny

    h = res.h
    dx = grid.dx
    dy = grid.dy
    ϕ = res.phi
    cw = fluid.cw
    co = fluid.co
    cf = res.cf
    dt = sim.dt
    s = sim.s
    pc = sim.pc

    D = zeros(Nx*Ny, Nx*Ny)
    d11 = zeros(Nx*Ny, Nx*Ny)
    d12 = zeros(Nx*Ny, Nx*Ny)
    d21 = zeros(Nx*Ny, Nx*Ny)
    d22 = zeros(Nx*Ny, Nx*Ny)

    for l in 1:Nx*Ny
        # j = ceil(Int64, l/Nx)
        # i = l - (j-l) * Nx
        V = h * dx * dy
        
        d11[l, l] = V * s[l] * ϕ * (cf + cw) / dt
        d12[l, l] = V * ϕ  / dt * (1 - s[l] * ϕ * cw * pc[l])
        d21[l,l] = V * ϕ * (1 - s[l]) / dt * (cf + co)
        d22[l,l] = - V / dt * ϕ

        D[l,l] = -d22[l,l] * d11[l,l] / d12[l,l] + d21[l,l]
    end

    return D, d11, d12, d21, d22
end