module IMPES

# Write your package code here.
include("types.jl")
include("flowFunctions.jl")
include("T_matrix.jl")
include("D_matrix.jl")
include("Q_vector.jl")


export RelPerms, Reservoir, Fluids, Grid, BoundaryConditions,
       Wells, Simulation

export D_matrix, T_matrix, Q_matrix
export IMPES!

end
