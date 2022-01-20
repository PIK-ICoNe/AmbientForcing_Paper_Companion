using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using CSV, DataFrames
using CairoMakie
using LightGraphs, FileIO
using BS_DAE
default(grid = false, foreground_color_legend = nothing, bar_edges = false, framestyle =:box, msc = :auto, dpi=300, legendfontsize = 11, labelfontsize = 12, tickfontsize = 10)

# Graph Plots and Trajectories in the Appendix
input = joinpath(@__DIR__,"../data/Pg_ensemble_method_Voltage_Networks_100_Nodes_100_SampleSize_250_TauMax1.0_Dist[[-2.0, 0.0], [-3.0, 2.0]].csv")
df = DataFrame(CSV.File(input))

# Voltage Drop Plot without desynronization
a = -428725142502246841
pg, rpg, op = find_power_grid(a)

plot_fault_simulation(pg, rpg, op, 50.0, [-1, 0, 0], 94)
plot_fault_simulation(pg, rpg, op, 50.0, [-1, 0, 0], 77)
plot_fault_simulation(pg, rpg, op, 50.0, [-2, -2, 0], 13)

f = my_graph_plot(pg, df, a, [94, 13, 77])

# Desynronization of 1 node and short circuit
c = 6845581020597587365
pg, rpg, op = find_power_grid(c)

plot_fault_simulation(pg, rpg, op, 1000.0, [-1, -1, 0], 62)

df[3262, :NodeClass]
f = my_graph_plot(pg, df, c, 62)

# Desynronization of 2 nodes and Voltage Drops
kappa = -8195088575577026753

pg, rpg, op = find_power_grid(kappa)
plot_fault_simulation(pg, rpg, op, 100.0, [-1, -0.5, 0], 45)

f = my_graph_plot(pg, df, kappa, 45)
df[2045, :NodeClass]