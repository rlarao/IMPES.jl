module IMPES

using LinearAlgebra, SparseArrays

include("types.jl")
include("flowFunctions.jl")
include("solvers.jl")
include("utilities.jl")
include("T_matrix.jl")

export RelPerms, SingleRelPerms, MixedRelPerms, Reservoir, Fluids, Grid, BoundaryConditions
export Simulation, FlowResults
export IMPES!
export CompTransport, IMPEC!
export recoveryFactor, injPoreVolumes, T_matrix
end


