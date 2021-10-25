function Thalf(l::Int64, m::Int64, res::Reservoir, kr::RelPerms,
            fluid::Fluids, grid::Grid, sim::Simulation)

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
    
    #Get i and j indices for blocks l and m
    lj= l ÷ Nx + 1
    li= l - (lj - 1) * Nx
    
    mj= l ÷ Nx + 1
    mi= m - (mj - 1) * Nx

    # Upwinding
    if p[l] >=p[m]
        krw, kro = kr(s[l], θ)
    elseif p[l] < p[m]
        krw, kro = kr(s[m], θ)
    end

    if lj==mj && li==mi # same cell
        if l-1 % Nx == 0 || l % Nx == 0
            Tw = k * Δy * h / (μw * Δx) * krw
            To = k * Δy * h / (μo * Δx) * kro
        else
            Tw = k * Δx * h / (μw * Δy) * krw
            To = k * Δx * h / (μo * Δy) * kro
        end
    
    elseif lj == mj # Get x-direction transmissibilities
        Tw = -  k * Δy * h / (μw * Δx) * krw
        To = -  k * Δy * h / (μo * Δx) * kro
    
    elseif li==mi # Get y-direction transmissibilities
        Tw = -  k * Δx * h / (μw * Δy) * krw
        To = -  k * Δx * h / (μo * Δy) * kro
    end

    return Tw, To
end


function T_matrix(res::Reservoir, kr::RelPerms, fluid::Fluids, grid::Grid,
                    sim::Simulation, bc::BoundaryConditions)
    
    Nx = grid.Nx
    Ny = grid.Ny
    
    # Initialize transmissibilities matrices
    To = zeros(Nx*Ny, Nx*Ny)
    Tw = zeros(Nx*Ny, Nx*Ny)

    for l in 1:Nx*Ny
        #* Left
        if (l-1) % Nx != 0 # block L is not on the left edge
            Tw[l, l-1], To[l, l-1] = Thalf(l, l-1, res, kr, fluid, grid, sim)
            Tw[l,l] -= Tw[l,l-1]
            To[l,l] -= To[l,l-1]
        elseif l-1 % Nx == 0 && bc.type[1] == :dirichlet
            Twb, Tob = Thalf(l,l, res, kr, fluid, grid, sim)
            Tw[l,l] -= 2*Twb;
            To[l,l] -= 2*Tob;
        end

        #* Right
        if l % Nx != 0 # block L is not on the right edge
            Tw[l, l+1], To[l, l+1] = Thalf(l, l+1, res, kr, fluid, grid, sim)
            Tw[l,l] -= Tw[l,l+1]
            To[l,l] -= To[l,l+1]
        elseif l % Nx == 0 && bc.type[2] == :dirichlet # block l is on the right edge with dirichlet BC 
            Twb, Tob =Thalf(l,l, res, fluid, grid, sim)
            Tw[l,l] -= 2*Twb;
            To[l,l] -= 2*Tob;
        end

        #* Bottom
        if ceil(Int64, l / Nx) != 1 # block L is not on the bottom edge
            Tw[l, l-Nx], To[l, l-Nx] = Thalf(l, l-Nx, res, kr,  fluid, grid, sim)
            Tw[l,l] -= Tw[l,l-Nx]
            To[l,l] -= To[l,l-Nx]
        elseif l % Nx == 0 && bc.type[3] == :dirichlet # block l is on the bottom edge with dirichlet BC 
            Twb, Tob =Thalf(l,l, res, kr, fluid, grid, sim)
            Tw[l,l] -= 2*Twb;
            To[l,l] -= 2*Tob;
        end

        #* Top
        if ceil(Int64, l / Nx) != Ny # block L is not on the bottom edge
            Tw[l, l+Nx], To[l, l+Nx] = Thalf(l, l+Nx, res, kr, fluid, grid, sim)
            Tw[l,l] -= Tw[l,l+Nx]
            To[l,l] -= To[l,l+Nx]
        elseif l % Nx == 0 && bc.type[4] == :dirichlet # block l is on the right edge with dirichlet BC 
            Twb, Tob =Thalf(l,l, res, kr, fluid, grid, sim)
            Tw[l,l] -= 2*Twb;
            To[l,l] -= 2*Tob;
        end
    end

    # T = -(sim.d12 \ sim.d22)* Tw + To
    T = Tw + To

    return T, Tw 
end

