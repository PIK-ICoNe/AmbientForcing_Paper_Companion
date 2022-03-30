using Pkg, Distributed
Pkg.activate(joinpath(@__DIR__, "../../"))
Pkg.instantiate()

using SlurmClusterManager
addprocs(SlurmManager())

@everywhere using Distributed, Pkg
@everywhere Pkg.activate(joinpath(@__DIR__, "../../"))
@everywhere using SlurmClusterManager

# All workers need to know the functions for the Computation
using BS_DAE
@everywhere using BS_DAE

input = joinpath(@__DIR__,"../../data/Pg_ensemble_method_Voltage_Networks_100_Nodes_100_SampleSize_250_TauMax1.0_Dist[[-2.0, 0.0], [-3.0, 2.0]].csv")
df = DataFrame(CSV.File(input))

@everywhere input = joinpath(@__DIR__,"../../data/Pg_ensemble_method_Voltage_Networks_100_Nodes_100_SampleSize_250_TauMax1.0_Dist[[-2.0, 0.0], [-3.0, 2.0]].csv")
@everywhere df = DataFrame(CSV.File(input))


# Node with most frequent desynchronization and no voltage drop
_, node_idx = findmax(df[!, :SNBS_u] .- df[!, :SNBS_ω]) 
pg_idx = df[node_idx, :PG_IDX] # related power grid
@everywhere _, node_idx = findmax(df[!, :SNBS_u] .- df[!, :SNBS_ω]) 
@everywhere pg_idx = df[node_idx, :PG_IDX] # related power grid

pg, rpg, op = find_power_grid(pg_idx)
@everywhere pg, rpg, op = find_power_grid(pg_idx)

parallel_snbs_surv(pg, op, 500.0, 100, [[-2.0, 0], [-3.0, 2.0]], 0.1, 1.0, "Voltage", false; write_frand = true)
