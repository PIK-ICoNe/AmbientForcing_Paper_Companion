## Compares different simulation times and solver tolarances!
using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using CSV, DataFrames
using BS_DAE

input_400 = joinpath(@__DIR__,"../data/Pg_ensemble_method_Voltage_Networks_100_Nodes_100_SampleSize_250_TauMax1.0_Dist[[-2.0, 0.0], [-3.0, 2.0]]_simtime_400.0.csv")
input_500 = joinpath(@__DIR__,"../data/Pg_ensemble_method_Voltage_Networks_100_Nodes_100_SampleSize_250_TauMax1.0_Dist[[-2.0, 0.0], [-3.0, 2.0]].csv")

df = DataFrame(CSV.File(input_400))
df_500 = DataFrame(CSV.File(input_500))

df_cons_500 = filter("Constraint" => x -> x == 1, df_500)
df_dyn_500 = filter("Constraint" => x -> x == 0, df_500)

df_cons = filter("Constraint" => x -> x == 1, df)
df_dyn = filter("Constraint" => x -> x == 0, df)

##########################################
#            Density Plots               #
##########################################

StatsPlots.density(df[!, :SNBS], lw = 3, label = "400s", legend = :topleft)
StatsPlots.density!(df_500[!, :SNBS], lw = 3, label = "500s")

StatsPlots.density(df[!, :SURV],  lw = 3, label = "400s", legend = :topleft)
StatsPlots.density!(df_500[!, :SURV], lw = 3, label = "500s")

################################################
# Voltage Drops, Infeasible, desynronization   #
################################################

p_vol, p_desyn, p_vol_only, p_desyn_only, p_infeas = snbs_fault_statistics(df)
snbs_fault_statistics(df_500)

########################################################
#                   Degree Plots                       #
########################################################
using StatsBase

SNBS_d_cons, SNBS_δ_cons, SURV_d_cons, SURV_δ_cons = filter_degree(df_cons)
SNBS_d_dyn, SNBS_δ_dyn, SURV_d_dyn, SURV_δ_dyn = filter_degree(df_dyn)

SNBS_d_cons_500, SNBS_δ_cons_500, SURV_d_cons_500, SURV_δ_cons_500 = filter_degree(df_cons_500)
SNBS_d_dyn_500, SNBS_δ_dyn_500, SURV_d_dyn_500, SURV_δ_dyn_500 = filter_degree(df_dyn_500)


Plots.plot(string.(collect(1:length(SNBS_d_cons))), SNBS_d_cons, ribbon = (SNBS_δ_cons[:,1], SNBS_δ_cons[:,2]), fillalpha = 0.3, lw = 3, color = ColorSchemes.leonardo[0.5], xlabel = L"d", ylabel = L"\beta", label = "Algebraic")
Plots.plot!(string.(collect(1:length(SNBS_d_dyn))), SNBS_d_dyn, ribbon = (SNBS_δ_dyn[:,1], SNBS_δ_dyn[:,2]), fillalpha = 0.3, lw = 3, color = ColorSchemes.leonardo[end], label = "Differential", legend = :bottomright)
Plots.plot!(string.(collect(1:length(SNBS_d_cons_500))), SNBS_d_cons_500, ribbon = (SNBS_δ_cons_500[:,1], SNBS_δ_cons_500[:,2]), fillalpha = 0.3, lw = 3, color = ColorSchemes.leonardo[0.5], xlabel = L"d", ylabel = L"\beta", label = "Algebraic 500")
Plots.plot!(string.(collect(1:length(SNBS_d_dyn_500))), SNBS_d_dyn_500, ribbon = (SNBS_δ_dyn_500[:,1], SNBS_δ_dyn_500[:,2]), fillalpha = 0.3, lw = 3, color = ColorSchemes.leonardo[end], label = "Differential 500", legend = :bottomright)

Plots.plot(string.(collect(1:length(SURV_d_cons))), lw = 3, SURV_d_cons, ribbon = (SURV_δ_cons[:,1], SURV_δ_cons[:,2]), color = ColorSchemes.berlin[0.25], xlabel = L"d", ylabel = L"\sigma", label = "Algebraic")
Plots.plot!(string.(collect(1:length(SURV_d_dyn))), lw = 3, SURV_d_dyn, ribbon = (SURV_δ_dyn[:,1], SURV_δ_dyn[:,2]), color = ColorSchemes.berlin[0.95], label = "Differential", legend = :bottomright)
Plots.plot!(string.(collect(1:length(SURV_d_cons_500))), lw = 3, SURV_d_cons_500, ribbon = (SURV_δ_cons_500[:,1], SURV_δ_cons_500[:,2]), color = ColorSchemes.berlin[0.25], xlabel = L"d", ylabel = L"\sigma", label = "Algebraic")
Plots.plot!(string.(collect(1:length(SURV_d_dyn_500))), lw = 3, SURV_d_dyn_500, ribbon = (SURV_δ_dyn_500[:,1], SURV_δ_dyn_500[:,2]), color = ColorSchemes.berlin[0.95], label = "Differential", legend = :bottomright)


density_degree_plot(df_dyn, "SURV", L"\sigma")
density_degree_plot!(df_dyn_500, "SURV", L"\sigma")

density_degree_plot(df_cons, "SURV", L"\sigma", xlims = [0.3, 1])
density_degree_plot!(df_cons_500, "SURV", L"\sigma", xlims = [0.3, 1])

density_degree_plot(df_dyn, "SNBS", L"\beta", color = ColorSchemes.leonardo, xlims = [0.5, 1])
density_degree_plot!(df_dyn_500, "SNBS", L"\beta", color = ColorSchemes.leonardo, xlims = [0.5, 1])


density_degree_plot(df_cons, "SNBS", L"\beta", color = ColorSchemes.leonardo, xlims = [0.75, 1])
density_degree_plot!(df_cons_500, "SNBS", L"\beta", color = ColorSchemes.leonardo, xlims = [0.75, 1])

# Tolerance



input_tol = joinpath(@__DIR__,"../data/Pg_ensemble_method_Voltage_Networks_100_Nodes_100_SampleSize_250_TauMax1.0_Dist[[-2.0, 0.0], [-3.0, 2.0]]_tol_1.0e-6.csv")

df = DataFrame(CSV.File(input_tol))
df_cons = filter("Constraint" => x -> x == 1, df)
df_dyn = filter("Constraint" => x -> x == 0, df)

##########################################
#            Density Plots               #
##########################################

StatsPlots.density(df[!, :SNBS], lw = 3, label = "Default Tol", legend = :topleft)
StatsPlots.density!(df_500[!, :SNBS], lw = 3, label = "1e-6")

StatsPlots.density(df[!, :SURV],  lw = 3, label = "Default Tol", legend = :topleft)
StatsPlots.density!(df_500[!, :SURV], lw = 3, label = "1e-6")


################################################
# Voltage Drops, Infeasible, desynronization   #
################################################

snbs_fault_statistics(df)
snbs_fault_statistics(df_500)

########################################################
#                   Degree Plots                       #
########################################################
using StatsBase

SNBS_d_cons, SNBS_δ_cons, SURV_d_cons, SURV_δ_cons = filter_degree(df_cons)
SNBS_d_dyn, SNBS_δ_dyn, SURV_d_dyn, SURV_δ_dyn = filter_degree(df_dyn)

Plots.plot(string.(collect(1:length(SNBS_d_cons))), SNBS_d_cons, ribbon = (SNBS_δ_cons[:,1], SNBS_δ_cons[:,2]), fillalpha = 0.3, lw = 3, color = ColorSchemes.leonardo[0.5], xlabel = L"d", ylabel = L"\beta", label = "Algebraic")
Plots.plot!(string.(collect(1:length(SNBS_d_dyn))), SNBS_d_dyn, ribbon = (SNBS_δ_dyn[:,1], SNBS_δ_dyn[:,2]), fillalpha = 0.3, lw = 3, color = ColorSchemes.leonardo[end], label = "Differential", legend = :bottomright)
Plots.plot!(string.(collect(1:length(SNBS_d_cons_500))), SNBS_d_cons_500, ribbon = (SNBS_δ_cons_500[:,1], SNBS_δ_cons_500[:,2]), fillalpha = 0.3, lw = 3, color = ColorSchemes.leonardo[0.5], xlabel = L"d", ylabel = L"\beta", label = "Algebraic 500")
Plots.plot!(string.(collect(1:length(SNBS_d_dyn_500))), SNBS_d_dyn_500, ribbon = (SNBS_δ_dyn_500[:,1], SNBS_δ_dyn_500[:,2]), fillalpha = 0.3, lw = 3, color = ColorSchemes.leonardo[end], label = "Differential 500", legend = :bottomright)



Plots.plot(string.(collect(1:length(SURV_d_cons))), lw = 3, SURV_d_cons, ribbon = (SURV_δ_cons[:,1], SURV_δ_cons[:,2]), color = ColorSchemes.berlin[0.25], xlabel = L"d", ylabel = L"\sigma", label = "Algebraic")
Plots.plot!(string.(collect(1:length(SURV_d_dyn))), lw = 3, SURV_d_dyn, ribbon = (SURV_δ_dyn[:,1], SURV_δ_dyn[:,2]), color = ColorSchemes.berlin[0.95], label = "Differential", legend = :bottomright)
Plots.plot!(string.(collect(1:length(SURV_d_cons_500))), lw = 3, SURV_d_cons_500, ribbon = (SURV_δ_cons_500[:,1], SURV_δ_cons_500[:,2]), color = ColorSchemes.berlin[0.25], xlabel = L"d", ylabel = L"\sigma", label = "Algebraic")
Plots.plot!(string.(collect(1:length(SURV_d_dyn_500))), lw = 3, SURV_d_dyn_500, ribbon = (SURV_δ_dyn_500[:,1], SURV_δ_dyn_500[:,2]), color = ColorSchemes.berlin[0.95], label = "Differential", legend = :bottomright)


density_degree_plot(df_dyn, "SURV", L"\sigma")
density_degree_plot!(df_dyn_500, "SURV", L"\sigma")

density_degree_plot(df_cons, "SURV", L"\sigma", xlims = [0.3, 1])
density_degree_plot!(df_cons_500, "SURV", L"\sigma", xlims = [0.3, 1])

density_degree_plot(df_dyn, "SNBS", L"\beta", color = ColorSchemes.leonardo, xlims = [0.5, 1])
density_degree_plot!(df_dyn_500, "SNBS", L"\beta", color = ColorSchemes.leonardo, xlims = [0.5, 1])


density_degree_plot(df_cons, "SNBS", L"\beta", color = ColorSchemes.leonardo, xlims = [0.75, 1])
density_degree_plot!(df_cons_500, "SNBS", L"\beta", color = ColorSchemes.leonardo, xlims = [0.75, 1])
