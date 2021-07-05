module IMPES

using LinearAlgebra

include("types.jl")
include("flowFunctions.jl")
include("T_matrix.jl")
include("D_matrix.jl")
include("Q_vector.jl")
include("solvers.jl")

export RelPerms, Reservoir, Fluids, Grid, BoundaryConditions,
       Wells, Simulation, PostProcessing
export D_matrix, T_matrix, Q_vector
export IMPES!

end
