using Pkg, Distributed
Pkg.activate(joinpath(@__DIR__, "../../"))
Pkg.instantiate()

using SlurmClusterManager, BS_DAE
addprocs(SlurmManager())

@everywhere using Distributed, Pkg
@everywhere Pkg.activate(joinpath(@__DIR__, "../../"))
@everywhere using SlurmClusterManager, BS_DAE

@everywhere num_nodes = 100
@everywhere num_networks = 100
@everywhere samples = 250

ensemble_calc([[-2, 0], [-3.0, 2.0]], 1.0, num_networks, samples,  num_nodes, "Voltage", spread_run = false, simtime = 400.0)