using Pkg, Distributed
Pkg.activate(joinpath(@__DIR__, "../../"))
Pkg.instantiate()
using SlurmClusterManager
addprocs(SlurmManager())

@everywhere using Distributed, Pkg
Pkg.activate(joinpath(@__DIR__, "../../"))
@everywhere using SlurmClusterManager

using BS_DAE
@everywhere using BS_DAE

@everywhere num_nodes = 100
@everywhere num_networks = 100
@everywhere samples = 250

ensemble_calc([[-2, 2], [-5.0, 5.0]], 1.0, num_networks, samples,  num_nodes, "All", spread_run = true)
