using Pkg, Distributed
Pkg.activate(joinpath(@__DIR__, "../../"))
Pkg.instantiate()

using SlurmClusterManager
addprocs(SlurmManager())

@everywhere using Distributed, Pkg
@everywhere Pkg.activate(joinpath(@__DIR__, "../../"))
@everywhere using SlurmClusterManager

# All workers need to know the functions for the Computation
using BS_DAE, PowerDynamics
@everywhere using BS_DAE, PowerDynamics

pg_ieee96 = get_ieee_96()
@everywhere pg_ieee96 = get_ieee_96()

@everywhere τ_max = 1.0
@everywhere op = find_operationpoint(pg_ieee96)

parallel_snbs_surv(pg_ieee96, op, 750.0, 500, [[-3.0, 3.0], [-3, 2]], 0.1, τ_max, "Voltage", true)
