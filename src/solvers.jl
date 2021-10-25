function IMPES!(sim::Simulation, res::Reservoir, kr::RelPerms, fluid::Fluids, grid::Grid,
                bc::BoundaryConditions, post::FlowResults)

    # Define flow functions
    # krw(s) = water_rel_perm(s, res.kr)
    # kro(s) = oil_rel_perm(s, res.kr)
    
    # λw(s) = krw(s) / fluid.μw
    # λo(s) = kro(s) / fluid.μo
    # λ(s) = λo(s) + λw(s)
    
    # fw(s) = λw(s) ./ (λw(s) .+ λo(s))

    # Get data
    n = grid.Nx
    A = res.A
    k = res.k
    ϕ = res.phi
    Δx = grid.dx
    Δt = sim.dt
    s = sim.s
    θ = sim.θ

    # Initialize velocity vectors
    v = zeros(n+1)
    v[1] = bc.value[1] / A
    
    ∇v = zeros(n)  

    # Calculate all mobilities
    krw, kro = kr(s, θ)
    Λw = krw ./ fluid.μw
    Λo = kro ./ fluid.μo
    Λ = Λw + Λo


    # Create T Vector
    T = zeros(n,n)
    for i in 2:n-1
        T[i,i-1] = Λ[i-1]
        T[i,i+1] = Λ[i]
        T[i,i] = - (T[i, i-1] + T[i,i+1])
    end

    # Left - Neumann
    T[1,1] = -Λ[1]
    T[1,2] = Λ[1]
    # Right - Dirichlet
    T[n,n] = - 9*Λ[n] - 3*Λ[n-1]
    T[n,n-1] = 3*Λ[n-1] + Λ[n]
    
    ## Get Q vector
    Q = zeros(n)
    # Left - Neumann
    Q[1] = - bc.value[1] / A / k * Δx
    # Right - Dirichlet
    Q[n] = -8 * (Λw[n] + Λo[n]) * bc.value[2] 

    ## Solve for Pressure    
    p = T \ Q

    ## Solve for Saturation
    # Calculate water velocity
    v[2:end-1] = [Λw[l-1] * (p[l-1] - p[l]) for l in 2:n] / Δx * k
    v[end] = 2 * Λw[n]*(p[n]-bc.value[2]) / Δx * k
    # Calculate gradient of water velocity
    ∇v = [v[l-1] - v[l] for l in 2:n+1] / Δx
    
    s = s .+ ∇v / ϕ * Δt
    
    # Update Simulation
    sim.p = p
    sim.s = s
    sim.v = v
    sim.i += 1
    sim.t += sim.dt

    # Save Results
    post.p[:,sim.i] = p
    post.s[:,sim.i] = s
    post.θ[:,sim.i] = θ

end


function IMPEC!(ctrans::CompTransport, sim, res, grid, post)
    n = grid.Nx
    Δx = grid.dx
    ϕ = res.phi
    Δt = sim.dt
    v = sim.v
    s = sim.s
    i = sim.i
    Δs = sim.s .- post.s[:,i-1]

    ncomps = ctrans.ncomps

    for j in 1:ncomps
        c = ctrans.c[:,j,i-1]

        cinj = ctrans.cinj[j]
        ∇cuw  = zeros(n)
        ∇cuw[2:n] = [c[l] * v[l+1] - c[l-1] * v[l]  for l in 2:n] / Δx
        ∇cuw[1] = (c[1] * v[2] - cinj * v[1])  / Δx

        c = c  .- Δt * ∇cuw / ϕ ./ s - c .* Δs ./ s

        ctrans.c[:,j,i] = c
    end

end



# function IMPEC!(ctrans::CompTransport, sim::Simulation,  res::Reservoir, fluid::Fluids, grid::Grid, bc::BoundaryConditions, wells::Wells)
#     Nx = grid.Nx
#     Ny = grid.Ny
#     p = sim.p
#     s = sim.s
#     ϕ = res.phi
#     k = res.k
#     Δx = grid.dx
#     Δt = sim.dt
#     μw = fluid.μw
#     Vi = res.h * res.W * Δx
#     qinj = wells.Q[1]
#     kr = res.kr
#     c = ctrans.c
#     cinj = ctrans.cinj
#     sold = ctrans.s_old

#     η = k * Δt / (ϕ * μw * Δx^2) * 6894757 * 1.06e-14 * 3600 * 24
#     A = zeros(Nx * Ny, Nx * Ny)
#     krw = [water_rel_perm(s[l], kr) for l in 1:Nx]

#     # Upwinding
#     for l in 2:Nx-1
#         # krw1 = water_rel_perm(s[l], kr)
#         # krw0 = water_rel_perm(s[l-1], kr)

#         A[l, l-1] =  krw[l-1] * (p[l-1] - p[l]) * η / s[l]
#         A[l, l] = krw[l] * (p[l+1] - p[l])  * η / s[l]
#     end
    
#     A[Nx, Nx-1] =  krw[Nx-1] * (p[Nx-1] - p[Nx]) * η / s[Nx]
#     A[Nx, Nx] =  - krw[Nx] * (2* p[Nx]) * η  / s[Nx]

#     A[1, 1] =  krw[1] * (p[2] - p[1]) * η  / s[1]
    
#     Q = zeros(Nx * Ny)
#     Q[1] = qinj * Δt / (s[1] * Vi * ϕ) * cinj
    
#     Δs = (s .- sold) ./ sold

#     ctrans.c = c + A * c + Q - Δs .* c
#     ctrans.s_old = sim.s
# end



# function IMPEC2!(ctrans::CompTransport, sim::Simulation,  res::Reservoir, fluid::Fluids, grid::Grid, bc::BoundaryConditions, wells::Wells)
#     Nx = grid.Nx
#     Ny = grid.Ny
#     p = sim.p
#     s = sim.s
#     ϕ = res.phi
#     k = res.k
#     Δx = grid.dx
#     Δt = sim.dt
#     μw = fluid.μw
#     Vi = res.h * res.W * Δx
#     qinj = wells.Q[1]
#     kr = res.kr
#     c = ctrans.c
#     cinj = ctrans.cinj
#     sold = ctrans.s_old

#     η = k * Δt / (ϕ * μw * Δx^2) * 6894757 * 1.06e-14 * 3600 * 24
#     A = zeros(Nx * Ny, Nx * Ny)

#     krw = [water_rel_perm(s[l], kr) for l in 1:Nx]
    
#     u = zeros(Nx)
#     for l in 2:Nx-1
#         u[l] = krw[l] * (p[l] - p[l+1])
#     end
#     u[1] = krw[1] * (p[2] - p[1])
#     u[Nx] = krw[Nx] * (p[Nx] - p[Nx-1])

#     # Upwinding
#     for l in 3:Nx - 1
#         A[l, l] = -3 * u[l] * η / 2 / s[l]
#         A[l, l-1] =  4 * u[l-1] * η / 2 / s[l]
#         A[l, l-2] =  -u[l-2] * η / 2 / s[l]
#     end
    
#     A[2, 1] = krw[1] * (p[1] - p[2]) * η / s[2]
#     A[2, 2] = krw[2] * (p[3] - p[2])  * η / s[2]

#     A[1, 1] =  krw[1] * (p[2] - p[1]) * η  / s[1]
    
#     A[Nx, Nx-1] =  krw[Nx-1] * (p[Nx-1] - p[Nx]) * η / s[Nx]
#     A[Nx, Nx] =  - krw[Nx] * (2* p[Nx]) * η  / s[Nx]

#     Q = zeros(Nx * Ny)
#     Q[1] = qinj * Δt / (s[1] * Vi * ϕ) * cinj
    
#     Δs = (s .- sold) ./ sold

#     ctrans.s_old = sim.s
#     ctrans.c = c + A * c + Q - Δs .* c
# end