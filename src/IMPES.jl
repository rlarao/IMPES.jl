module IMPES

using LinearAlgebra, SparseArrays

include("types.jl")
include("flowFunctions.jl")
include("solvers.jl")
# include("T_matrix.jl")
# include("D_matrix.jl")
# include("Q_vector.jl")

export RelPerms, Reservoir, Fluids, Grid, BoundaryConditions
export water_rel_perm, oil_rel_perm, Simulation, FlowResults# export D_matrix, T_matrix, Q_vector
export IMPES!
export CompTransport, IMPEC!

end


