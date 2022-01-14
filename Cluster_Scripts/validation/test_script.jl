using Pkg, Distributed
Pkg.activate(joinpath(@__DIR__, "../"))
Pkg.instantiate()
#addprocs()

#@everywhere using Distributed, Pkg
#@everywhere Pkg.activate(joinpath(@__DIR__, "../"))

# All workers need to know the functions for the Computation
include(joinpath(@__DIR__,"../../src/SNBS.jl"))
#@everywhere include(joinpath(@__DIR__,"../src/SNBS.jl"))

include(joinpath(@__DIR__,"../../src/helper_functions.jl"))
#@everywhere include(joinpath(@__DIR__,"../src/helper_functions.jl"))

#include(joinpath(@__DIR__,"../powergrids/ieee-rts-96-master/PowerDynamics_dyn.jl"))
#@everywhere include(joinpath(@__DIR__,"../powergrids/ieee-rts-96-master/PowerDynamics_dyn.jl"))

#@everywhere τ_max = 1.0
#@everywhere op = find_operationpoint(pg_ieee96)

#@everywhere dist_vec = [[-0.1, 0.1], [-3.0, 2.0]]
#parallel_snbs_surv(pg_ieee96, op, 1.0, 1, dist_vec, 0.1, τ_max, "All", false)


include(joinpath(@__DIR__, "../../src/pg_ensemble.jl"))
#@everywhere include(joinpath(@__DIR__, "../src/pg_ensemble.jl"))

#@everywhere num_nodes = 10
#@everywhere num_networks = 1
#@everywhere samples = 1


ensemble_calc([[-100, 100], [-3.0, 2.0]], 1.0, 1, 1, 10, "Voltage", spread_run = false, tol = 1e-6)