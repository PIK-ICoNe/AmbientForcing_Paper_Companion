using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using CSV, DataFrames, Plots
using BS_DAE, PowerDynamics
default(grid = false, foreground_color_legend = nothing, bar_edges = false, framestyle =:box, msc = :auto, dpi=300, legendfontsize = 11, labelfontsize = 12, tickfontsize = 10)

# Checking if nodes desynchronize without a voltage drop
input = joinpath(@__DIR__,"../data/Pg_ensemble_method_Voltage_Networks_100_Nodes_100_SampleSize_250_TauMax1.0_Dist[[-2.0, 0.0], [-3.0, 2.0]].csv")
df_ensemble = DataFrame(CSV.File(input))

_, node_idx = findmax(df_ensemble[!, :SNBS_u] .- df_ensemble[!, :SNBS_ω])
pg_idx = df_ensemble[node_idx, :PG_IDX] # related power grid
pg, rpg, op = find_power_grid(pg_idx)

ur_idx = findall(map(variable -> occursin("u_r", variable), string.(rpg.syms)))
ui_idx = findall(map(variable -> occursin("u_i", variable), string.(rpg.syms)))

input = joinpath(@__DIR__,"../data/Frands_Nodes_100.csv")
df = DataFrame(CSV.File(input))

# Node with most frequent desynchronization and no voltage drop
sim_time = 500.0

z_new = df[44, :F_RAND]
z_new = replace(z_new, "[" => "" )
z_new = replace(z_new, "]" => "" )
z_new = split.(z_new,",")
z_new = parse.(Float64, z_new)

sol = simulate_perturbation(rpg, z_new, (0.0, sim_time))
# We can see that there is in fact a deviation of the voltage from the fixed point
# This is also the case for the other networks
# We can confirm that the desynchronization never occurs without deviations of the voltage
plot_res(PowerGridSolution(sol, pg))

final_state = sol.u[end]
final_diff_v = abs.(op[:, :v] .- (abs.(final_state[ur_idx] .+ 1im * final_state[ui_idx]) ./ op[:, :v]))
maximum(final_diff_v)