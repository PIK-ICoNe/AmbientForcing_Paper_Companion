using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using CSV, DataFrames
using BS_DAE, Plots, Plots.Measures
using LaTeXStrings
using StatsBase, ColorSchemes
default(grid = false, foreground_color_legend = nothing, bar_edges = false, framestyle =:box, msc = :auto, dpi=300, legendfontsize = 11, labelfontsize = 12, tickfontsize = 10)

input = joinpath(@__DIR__,"../data/Pg_ensemble_method_Voltage_Networks_100_Nodes_100_SampleSize_250_TauMax1.0_Dist[[-2.0, 0.0], [-3.0, 2.0]].csv")
df = DataFrame(CSV.File(input))

df_cons = filter("Constraint" => x -> x == 1, df)
df_dyn = filter("Constraint" => x -> x == 0, df)

################################################
# Voltage Drops, Infeasible, desynronization   #
################################################
snbs_fault_statistics(df)

######################################
#     Topological Properties         #
######################################

density_class_plot(df_cons, "SNBS", L"\beta")
density_class_plot(df_dyn, "SNBS", L"\beta")
density_class_plot(df_cons, "SURV", L"\sigma")
density_class_plot(df_dyn, "SURV", L"\sigma", legend_pos = :topright)

########################################################
#                   Degree Plots                       #
########################################################

SNBS_d_cons, SNBS_δ_cons, SURV_d_cons, SURV_δ_cons = filter_degree(df_cons)
SNBS_d_dyn, SNBS_δ_dyn, SURV_d_dyn, SURV_δ_dyn = filter_degree(df_dyn)

Plots.plot(string.(collect(1:length(SNBS_d_cons))), SNBS_d_cons, ribbon = (SNBS_δ_cons[:,1], SNBS_δ_cons[:,2]), fillalpha = 0.3, lw = 3, color = ColorSchemes.leonardo[0.5], xlabel = L"d", ylabel = L"\beta", label = "Algebraic")
Plots.plot!(string.(collect(1:length(SNBS_d_dyn))), SNBS_d_dyn, ribbon = (SNBS_δ_dyn[:,1], SNBS_δ_dyn[:,2]), fillalpha = 0.3, lw = 3, color = ColorSchemes.leonardo[end], label = "Differential", legend = :bottomright)

Plots.plot(string.(collect(1:length(SURV_d_cons))), lw = 3, SURV_d_cons, ribbon = (SURV_δ_cons[:,1], SURV_δ_cons[:,2]), color = ColorSchemes.berlin[0.25], xlabel = L"d", ylabel = L"\sigma", label = "Algebraic")
Plots.plot!(string.(collect(1:length(SURV_d_dyn))), lw = 3, SURV_d_dyn, ribbon = (SURV_δ_dyn[:,1], SURV_δ_dyn[:,2]), color = ColorSchemes.berlin[0.95], label = "Differential", legend = :bottomright)

density_degree_plot(df_dyn, "SURV", L"\sigma", L"p(\sigma | d)", xlims = [0, 1])
density_degree_plot(df_cons, "SURV", L"\sigma", L"p(\sigma | d)", xlims = [0.3, 1])

#########################################
#           Infeasabilities             #
#########################################

jensen_shannon_divergence(df_cons[!, :SNBS], 1 .- df_cons[!, :INFEASIBLE], 0:0.001:1)
jensen_shannon_divergence(df_dyn[!, :SNBS], 1 .- df_dyn[!, :INFEASIBLE], 0:0.001:1)
jensen_shannon_divergence(df_cons[!, :SURV], 1 .- df_cons[!, :INFEASIBLE], 0:0.001:1)
jensen_shannon_divergence(df_dyn[!, :SURV], 1 .- df_dyn[!, :INFEASIBLE], 0:0.001:1)

#####################################
#   Power Output, no difference     #
#####################################
df[!, :Power]
df_p = filter("Power" => x -> x != "Slack", df) # filter out all Slacks
df_p[!, :Power] = map(x -> parse(Float64, x), df_p[!, :Power])

df_producer = filter("Power" => x -> x > 0, df_p) 
df_consumer = filter("Power" => x -> x < 0, df_p) 

@df df_producer StatsPlots.density(:SNBS, lw = 3, group = :Constraint, xlims = [0, 1], line = :dash, fillalpha = 0.8, bins = 12, label = ["Differential" "Algebraic"], legend = :topleft)
@df df_consumer StatsPlots.density!(:SNBS, lw = 3, group = :Constraint, xlims = [0, 1], fillalpha = 0.8, bins = 12, label = ["" ""])

@df df_producer StatsPlots.density(:SURV, lw = 3, group = :Constraint, xlims = [0, 1], line = :dash, fillalpha = 0.8, bins = 12, label = ["Differential" "Algebraic"], legend = :topleft)
@df df_consumer StatsPlots.density!(:SURV, lw = 3, group = :Constraint, xlims = [0, 1], fillalpha = 0.8, bins = 12, legend = :topleft, label = ["" ""])