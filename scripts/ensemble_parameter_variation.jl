# Comparing different Parameters for SchifferApprox

using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using CSV, DataFrames
using BS_DAE

input = joinpath(@__DIR__,"../data/Pg_ensemble_method_All_Networks_100_Nodes_100_SampleSize_100_Dist[[-2.0, 2.0], [-5.0, 5.0]]_5.5_8.8_5.5_0.11.csv")
df = DataFrame(CSV.File(input))

df_cons = filter("Constraint" => x -> x == 1, df)
df_dyn = filter("Constraint" => x -> x == 0, df)

################################################
# Voltage Drops, Infeasible, desynronization   #
################################################

p_vol, p_desyn, p_vol_only, p_desyn_only, p_infeas = snbs_fault_statistics(df)



############################################
#    Frequency vs. Voltage Conditions      #
############################################

@df df scatter(:SNBS_u, :SNBS_ω, group = :Constraint, yaxis = L"\beta_{\omega}", xaxis = L"\beta_u", legend = :bottomright, markersize = 5, label = ["Differential" "Algebraic"])
Plots.plot!([0,1], [0,1], line = :dash, label = "", color = colorant"gray0")

@df df scatter(:SURV_u, :SURV_ω, group = :Constraint, yaxis = L"\sigma_{\omega}", xaxis = L"\sigma_u", legend = :bottomright, markersize = 5, label = ["Differential" "Algebraic"])
Plots.plot!([0,1], [0,1], line = :dash, label = "", color = colorant"gray0")

@df df groupedhist(:SNBS, group = :Constraint, xlims = [0, 1], bins = 12, legend = :topleft, fillalpha = 0.7, color = [colorant"lightseagreen" colorant"firebrick1" ],bar_position = :stack,   xaxis = L"\beta_{\omega, u}", bar_edges = false, linewidth = 0, label = ["Differential" "Algebraic"])
@df df groupedhist(:SURV, group = :Constraint, xlims = [0, 1], bins = 12, legend = :topleft, fillalpha = 0.7, color = [colorant"deeppink" colorant"green" ],bar_position = :stack,  xaxis = L"\sigma_{\omega, u}", bar_edges = false, linewidth = 0, label = ["Differential" "Algebraic"])

########################################################
#                   Degree Plots                       #
########################################################
using StatsBase

SNBS_d_cons, SNBS_δ_cons, SURV_d_cons, SURV_δ_cons = filter_degree(df_cons)
SNBS_d_dyn, SNBS_δ_dyn, SURV_d_dyn, SURV_δ_dyn = filter_degree(df_dyn)

Plots.plot(string.(collect(1:length(SNBS_d_cons))), SNBS_d_cons, ribbon = (SNBS_δ_cons[:,1], SNBS_δ_cons[:,2]), fillalpha = 0.3, lw = 3, color = colorant"firebrick", xlabel = L"d", ylabel = L"\beta", label = "Algebraic")
Plots.plot!(string.(collect(1:length(SNBS_d_dyn))), SNBS_d_dyn, ribbon = (SNBS_δ_dyn[:,1], SNBS_δ_dyn[:,2]), fillalpha = 0.3, lw = 3, color = colorant"coral1", label = "Differential", legend = :bottomright)

Plots.plot(string.(collect(1:length(SURV_d_cons))), lw = 3, SURV_d_cons, ribbon = (SURV_δ_cons[:,1], SURV_δ_cons[:,2]), color = colorant"teal", xlabel = L"d", ylabel = L"\sigma", label = "Algebraic")
Plots.plot!(string.(collect(1:length(SURV_d_dyn))), lw = 3, SURV_d_dyn, ribbon = (SURV_δ_dyn[:,1], SURV_δ_dyn[:,2]), color = colorant"steelblue", label = "Differential", legend = :bottomright)

#########################################
#           Infeasabilities             #
#########################################

StatsPlots.density(abs.(df_cons[!, :SNBS] .- (1 .-  df_cons[!, :INFEASIBLE])), lw = 3, label = "Algebraic")
StatsPlots.density!(abs.(df_dyn[!, :SNBS] .- (1 .- df_dyn[!, :INFEASIBLE])), lw = 3, label = "Differential", xaxis = L"\beta - f")

StatsPlots.density(abs.(df_cons[!, :SURV] .- (1 .-  df_cons[!, :INFEASIBLE])), lw = 3, label = "Algebraic")
StatsPlots.density!(abs.(df_dyn[!, :SURV] .- (1 .- df_dyn[!, :INFEASIBLE])), lw = 3, label = "Differential", xaxis = L"\sigma - f")


scatter(1 .- df_dyn[!, :INFEASIBLE], df_dyn[!, :SNBS], yaxis = L"\beta", xaxis = L"f", markersize = 5, label = "Algebraic")
scatter!(1 .- df_cons[!, :INFEASIBLE], df_cons[!, :SNBS], yaxis = L"\beta", xaxis = L"f", legend = :topleft, markersize = 5, label = "Differential")
Plots.plot!([0,1], [0,1], line = :dash, label = "", color = colorant"gray0")


scatter(1 .- df_cons[!, :INFEASIBLE], df_cons[!, :SURV],  yaxis = L"\sigma", xaxis = L"f", markersize = 5, label = "Algebraic")
scatter!(1 .- df_dyn[!, :INFEASIBLE], df_dyn[!, :SURV], legend = :topleft, yaxis = L"\sigma", xaxis = L"f", markersize = 5, label = "Differential")
