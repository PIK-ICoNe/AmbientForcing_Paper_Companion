using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using CSV, DataFrames
using CairoMakie
using GraphMakie
using LightGraphs, FileIO
using BS_DAE

input = joinpath(@__DIR__,"../data/Pg_ensemble_method_Voltage_Networks_100_Nodes_100_SampleSize_250_TauMax1.0_Dist[[-2.0, 0.0], [-3.0, 2.0]].csv")
df = DataFrame(CSV.File(input))

# Voltage Drop Plots

a = -428725142502246841
pg, rpg, op = find_power_grid(a)

plot_fault_simulation(pg, rpg, op, 50.0, [-1, 0, 0], 94)
png(joinpath(@__DIR__,"../plots/ensemble/traj/trajectory_node_94_pg_idx_" * string(a)))

plot_fault_simulation(pg, rpg, op, 50.0, [-1, 0, 0], 77)
png(joinpath(@__DIR__,"../plots/ensemble/traj/trajectory_node_77_pg_idx_" * string(a)))

plot_fault_simulation(pg, rpg, op, 50.0, [-2, -2, 0], 13)
png(joinpath(@__DIR__,"../plots/ensemble/traj/trajectory_node_13_pg_idx_" * string(a)))

f = my_graph_plot(pg, df, a, [94, 13, 77])
plot_place = joinpath(@__DIR__,"../plots/ensemble/traj/", "graph_snbs_pg_idx_" * string(a) * ".svg")
#FileIO.save(plot_place, f, resolution = (1000, 1000))

# Desynronization of 1 node
c = 6845581020597587365
pg, rpg, op = find_power_grid(c)

plot_fault_simulation(pg, rpg, op, 1000.0, [-1, -1, 0], 62)
png(joinpath(@__DIR__,"../plots/ensemble/traj/trajectory_node_62_pg_idx_" * string(c)))

df[3262, :NodeClass]
f = my_graph_plot(pg, df, c, 62)
plot_place = joinpath(@__DIR__,"../plots/ensemble/traj/", "graph_snbs_pg_idx_" * string(c) * ".svg")
#FileIO.save(plot_place, f, resolution = (1000, 1000))


# Desynronization of 2 nodes
kappa = -8195088575577026753

pg, rpg, op = find_power_grid(kappa)
plot_fault_simulation(pg, rpg, op, 100.0, [-1, -0.5, 0], 45)
png(joinpath(@__DIR__,"../plots/ensemble/traj/trajectory_node_45_pg_idx_" * string(kappa)))


f = my_graph_plot(pg, df, kappa, 45)
plot_place = joinpath(@__DIR__,"../plots/ensemble/traj/", "graph_snbs_pg_idx_" * string(kappa) * ".svg")
#FileIO.save(plot_place, f, resolution = (1000, 1000))

df[2045, :NodeClass]

#
b = -7880086776645138134
pg, rpg, op = find_power_grid(b)
plot_fault_simulation(pg, rpg, op, 100.0, [-0.5, -0.5, 0], 20)
png(joinpath(@__DIR__,"../plots/ensemble/traj/trajectory_node_20_pg_idx_" * string(b)))

f = my_graph_plot(pg, df, b, 20)
plot_place = joinpath(@__DIR__,"../plots/ensemble/traj/", "graph_snbs_pg_idx_" * string(b) * ".svg")
#FileIO.save(plot_place, f, resolution = (1000, 1000))