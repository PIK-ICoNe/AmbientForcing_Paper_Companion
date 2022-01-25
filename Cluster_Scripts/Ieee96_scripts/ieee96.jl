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

@everywhere dist_vec = [[-0.1, 0.1], [-3.0, 2.0]]
parallel_snbs_surv(pg_ieee96, op, 1000.0, 500, dist_vec, 0.1, τ_max, "All", true)

@everywhere dist_vec = [[-0.25, 0.25], [-5.0, 5.0]]
parallel_snbs_surv(pg_ieee96, op, 1000.0, 500, dist_vec, 0.1, τ_max, "All", true)

@everywhere dist_vec = [[-0.5, 0.5], [-10.0, 10.0]]
parallel_snbs_surv(pg_ieee96, op, 1000.0, 500, dist_vec, 0.1, τ_max, "All", true)