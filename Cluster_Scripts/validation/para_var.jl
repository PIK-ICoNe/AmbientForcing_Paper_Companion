using Pkg, Distributed
Pkg.activate(joinpath(@__DIR__, "../"))
Pkg.instantiate()
using SlurmClusterManager
addprocs(SlurmManager())

@everywhere using Distributed, Pkg
@everywhere Pkg.activate(joinpath(@__DIR__, "../"))
@everywhere using SlurmClusterManager

include(joinpath(@__DIR__, "../../src/pg_ensemble.jl"))
@everywhere include(joinpath(@__DIR__, "../../src/pg_ensemble.jl"))

@everywhere num_nodes = 100
@everywhere num_networks = 100
@everywhere samples = 100

ensemble_calc([[-2, 2], [-5.0, 5.0]], 1.0, num_networks, samples,  num_nodes, "All", τ_P = 5.5, τ_Q = 8.8, K_P = 5.5, K_Q = 0.11)
