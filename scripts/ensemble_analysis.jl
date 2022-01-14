using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using CSV, DataFrames
using BS_DAE

input = joinpath(@__DIR__,"../data/Pg_ensemble_method_All_Networks_100_Nodes_100_SampleSize_250_TauMax1.0_Dist[[-2.0, 2.0], [-5.0, 5.0]].csv")
df = DataFrame(CSV.File(input))

df_cons = filter("Constraint" => x -> x == 1, df)
df_dyn = filter("Constraint" => x -> x == 0, df)

################################################
# Voltage Drops, Infeasible, desynronization   #
################################################

p_vol, p_desyn, p_vol_only, p_desyn_only, p_infeas = snbs_fault_statistics(df)

######################################
#     Topological Properties         #
######################################

@df df Plots.scatter(:SNBS, :SURV, group = :NodeClass, xaxis = L"\beta", yaxis = L"\sigma", xlims = [0, 1], ylims = [0,1], legend = :topleft)
Plots.plot!([0,1], [0,1], line = :dash, label = "", color = colorant"gray0")
#png(joinpath(@__DIR__,"../plots/ensemble/scatter_topological"))

density_class_plot(df_cons, "SNBS", L"\beta")
#png(joinpath(@__DIR__,"../plots/ensemble/snbs_dens_topological_cons"))

density_class_plot(df_dyn, "SNBS", L"\beta")
#png(joinpath(@__DIR__,"../plots/ensemble/snbs_dens_topological_dyn"))

density_class_plot(df_cons, "SURV", L"\sigma")
#png(joinpath(@__DIR__,"../plots/ensemble/surv_dens_topological_cons"))

density_class_plot(df_dyn, "SURV", L"\sigma", legend_pos = :topright)
#png(joinpath(@__DIR__,"../plots/ensemble/surv_dens_topological_dyn"))

############################################
#    Frequency vs. Voltage Conditions      #
############################################

@df df Plots.scatter(:SNBS_u, :SNBS_ω, group = :Constraint, yaxis = L"\beta_{\omega}", xaxis = L"\beta_u", legend = :bottomright, markersize = 5, label = ["Differential" "Algebraic"])
Plots.plot!([0,1], [0,1], line = :dash, label = "", color = colorant"gray0")
#png(joinpath(@__DIR__,"../plots/ensemble/ieee96/snbs_u_omega"))

@df df Plots.scatter(:SURV_u, :SURV_ω, group = :Constraint, yaxis = L"\sigma_{\omega}", xaxis = L"\sigma_u", legend = :bottomright, markersize = 5, label = ["Differential" "Algebraic"])
Plots.plot!([0,1], [0,1], line = :dash, label = "", color = colorant"gray0")
#png(joinpath(@__DIR__,"../plots/ensemble/surv_u_omega"))


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
##png(joinpath(@__DIR__,"../plots/ensemble/snbs_degree"))

Plots.plot(string.(collect(1:length(SURV_d_cons))), lw = 3, SURV_d_cons, ribbon = (SURV_δ_cons[:,1], SURV_δ_cons[:,2]), color = colorant"teal", xlabel = L"d", ylabel = L"\sigma", label = "Algebraic")
Plots.plot!(string.(collect(1:length(SURV_d_dyn))), lw = 3, SURV_d_dyn, ribbon = (SURV_δ_dyn[:,1], SURV_δ_dyn[:,2]), color = colorant"steelblue", label = "Differential", legend = :bottomright)
##png(joinpath(@__DIR__,"../plots/ensemble/surv_degree"))

#########################################
#           Infeasabilities             #
#########################################
using KernelDensity
using StatsBase

density_estimator_cons = kde(abs.(df_cons[!, :SNBS] .- (1 .-  df_cons[!, :INFEASIBLE])))
max_cons, pos_cons = findmax(density_estimator_cons.density)

density_estimator_dyn = kde(abs.(df_dyn[!, :SNBS] .- (1 .-  df_dyn[!, :INFEASIBLE])))
max_dyn, pos_dyn = findmax(density_estimator_dyn.density)

plot(density_estimator_cons.x, density_estimator_cons.density / max_cons, lw = 3, label = "Algebraic")
plot!(density_estimator_dyn.x, density_estimator_dyn.density / max_dyn, lw = 3, label = "Differential", xaxis = L"\beta - f")
png(joinpath(@__DIR__,"../plots/ensemble/snbs_infeasible_dens"))


density_estimator_cons = kde(abs.(df_cons[!, :SURV] .- (1 .-  df_cons[!, :INFEASIBLE])))
max_cons, pos_cons = findmax(density_estimator_cons.density)

density_estimator_dyn = kde(abs.(df_dyn[!, :SURV] .- (1 .- df_dyn[!, :INFEASIBLE])))
max_dyn, pos_dyn = findmax(density_estimator_dyn.density)

plot(density_estimator_cons.x, density_estimator_cons.density / max_cons, lw = 3, label = "Algebraic")
plot!(density_estimator_dyn.x, density_estimator_dyn.density / max_dyn, lw = 3, label = "Differential", xaxis = L"\sigma - f")
png(joinpath(@__DIR__,"../plots/ensemble/surv_infeasible"))

########################################
#   Power FLOW,                        #
########################################

df[!, :Power]
df_p = filter("Power" => x -> x != "Slack", df) # filter out all Slacks
df_p[!, :Power] = map(x -> parse(Float64, x), df_p[!, :Power])


df_producer = filter("Power" => x -> x > 0, df_p) 
df_consumer = filter("Power" => x -> x < 0, df_p) 

@df df_producer StatsPlots.density(:SNBS, lw = 3, group = :Constraint, xlims = [0, 1], line = :dash, fillalpha = 0.8, bins = 12, label = ["Differential" "Algebraic"], legend = :topleft)
@df df_consumer StatsPlots.density!(:SNBS, lw = 3, group = :Constraint, xlims = [0, 1], fillalpha = 0.8, bins = 12, label = ["" ""])

@df df_producer StatsPlots.density(:SURV, lw = 3, group = :Constraint, xlims = [0, 1], line = :dash, fillalpha = 0.8, bins = 12, legend = false)
@df df_consumer StatsPlots.density!(:SURV, lw = 3, group = :Constraint, xlims = [0, 1], fillalpha = 0.8, bins = 12, legend = :topleft)


Plots.plot(df[!, "betweenness"] .* 100, df[!, "SNBS"], seriestype =:scatter, legend = false)
Plots.plot(df[!, "betweenness"] .* 100, df[!, "SURV"], seriestype =:scatter, legend = false)

#############################
#    Spreadability Plots    #
#############################

nan_idx = findall(map(x -> isnan(x), df[!, :Abs_Spread]))
delete!(df, nan_idx)

spread_by_dist, spread_sorted, distances_sorted = sort_by_distance(df[!, "Abs_Spread"], df[!, "Constraint_Distance"])
dist_zero =  df[!, "Constraint_Distance"] .* abs.(1 .- df[!, "Constraint"])
spread_by_dist_0, spread_sorted_0, distances_sorted_0 = sort_by_distance(df[!, "Abs_Spread"], dist_zero)

plt1 = spread_violin_plot(distances_sorted, spread_sorted, colorant"orange2", colorant"steelblue", latexstring(string(raw"\min (d", "_{cn})")))
plt2 = spread_violin_plot(distances_sorted_0, spread_sorted_0, colorant"peachpuff3", colorant"royalblue1", latexstring(string(raw"\min (d", "_{c})")))

savefig(plt1, joinpath(@__DIR__, "../plots/ensemble/spread_dcn.png"))
savefig(plt2, joinpath(@__DIR__, "../plots/ensemble/spread_dc.png"))

###############################
#  Amplification vs. Degree   #
###############################
using KernelDensity, StatsBase, Dierckx

Ax_degree, Ax_d_δ = filter_degree_Ax(df)
my_range = 0:0.01:4

dens = kde(df[!,:AX])
dens = Spline1D(dens.x, dens.density)

plot(my_range, dens(my_range) ./ sum(dens(my_range)),  lw = 3, xaxis = L"A", yaxis = L"p(A)", legend = false)
png(joinpath(@__DIR__, "../plots/ensemble/Amplification_dens"))

Plots.plot(string.(collect(1:length(Ax_degree))), Ax_degree, ribbon = (Ax_d_δ[:,1], Ax_d_δ[:,2]), fillalpha = 0.3, lw = 3, xaxis = L"d", yaxis = L"A", legend = false)
png(joinpath(@__DIR__, "../plots/ensemble/Amplification_degree"))
