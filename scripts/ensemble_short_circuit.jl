using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using CSV, DataFrames
using BS_DAE, Plots, Plots.Measures
using LaTeXStrings

input = joinpath(@__DIR__,"../data/Pg_ensemble_method_Voltage_Networks_100_Nodes_100_SampleSize_250_TauMax1.0_Dist[[-2.0, 0.0], [-3.0, 2.0]].csv")
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

density_class_plot(df_cons, "SNBS", L"\beta")
png(joinpath(@__DIR__,"../plots/ensemble/short_circ/snbs_dens_topological_cons"))

density_class_plot(df_dyn, "SNBS", L"\beta")
png(joinpath(@__DIR__,"../plots/ensemble/short_circ/snbs_dens_topological_dyn"))

density_class_plot(df_cons, "SURV", L"\sigma")
png(joinpath(@__DIR__,"../plots/ensemble/short_circ/surv_dens_topological_cons"))

density_class_plot(df_dyn, "SURV", L"\sigma", legend_pos = :topright)
png(joinpath(@__DIR__,"../plots/ensemble/short_circ/surv_dens_topological_dyn"))

########################################################
#                   Degree Plots                       #
########################################################
using StatsBase, ColorSchemes

SNBS_d_cons, SNBS_δ_cons, SURV_d_cons, SURV_δ_cons = filter_degree(df_cons)
SNBS_d_dyn, SNBS_δ_dyn, SURV_d_dyn, SURV_δ_dyn = filter_degree(df_dyn)

Plots.plot(string.(collect(1:length(SNBS_d_cons))), SNBS_d_cons, ribbon = (SNBS_δ_cons[:,1], SNBS_δ_cons[:,2]), fillalpha = 0.3, lw = 3, color = ColorSchemes.leonardo[0.5], xlabel = L"d", ylabel = L"\beta", label = "Algebraic")
Plots.plot!(string.(collect(1:length(SNBS_d_dyn))), SNBS_d_dyn, ribbon = (SNBS_δ_dyn[:,1], SNBS_δ_dyn[:,2]), fillalpha = 0.3, lw = 3, color = ColorSchemes.leonardo[end], label = "Differential", legend = :bottomright)
png(joinpath(@__DIR__,"../plots/ensemble/short_circ/snbs_degree"))

Plots.plot(string.(collect(1:length(SURV_d_cons))), lw = 3, SURV_d_cons, ribbon = (SURV_δ_cons[:,1], SURV_δ_cons[:,2]), color = ColorSchemes.berlin[0.25], xlabel = L"d", ylabel = L"\sigma", label = "Algebraic")
Plots.plot!(string.(collect(1:length(SURV_d_dyn))), lw = 3, SURV_d_dyn, ribbon = (SURV_δ_dyn[:,1], SURV_δ_dyn[:,2]), color = ColorSchemes.berlin[0.95], label = "Differential", legend = :bottomright)
png(joinpath(@__DIR__,"../plots/ensemble/short_circ/surv_degree"))

den = density_degree_plot(df_dyn, "SURV", L"\sigma", L"p(\sigma | d)", xlims = [0, 1])
png(joinpath(@__DIR__,"../plots/ensemble/short_circ/surv_dens_degree_dyn"))

density_degree_plot(df_cons, "SURV", L"\sigma", L"p(\sigma | d)", xlims = [0.3, 1])
png(joinpath(@__DIR__,"../plots/ensemble/short_circ/surv_dens_degree_cons"))

density_degree_plot(df_dyn, "SNBS", L"\beta", L"p(\sigma | d)",color = ColorSchemes.leonardo, xlims = [0.5, 1])
png(joinpath(@__DIR__,"../plots/ensemble/short_circ/snbs_dens_degree_dyn"))

density_degree_plot(df_cons, "SNBS", L"\beta", L"p(\sigma | d)", color = ColorSchemes.leonardo, xlims = [0.75, 1])
png(joinpath(@__DIR__,"../plots/ensemble/short_circ/snbs_dens_degree_cons"))

############################################
#    Frequency vs. Voltage Conditions      #
############################################

@df df Plots.scatter(:SNBS_u, :SNBS_ω, group = :Constraint, yaxis = L"\beta_{\omega}", xaxis = L"\beta_u", legend = :bottomright, markersize = 5, label = ["Differential" "Algebraic"])
Plots.plot!([0,1], [0,1], line = :dash, label = "", color = colorant"gray0")
png(joinpath(@__DIR__,"../plots/ensemble/short_circ/snbs_u_omega"))

@df df Plots.scatter(:SURV_u, :SURV_ω, group = :Constraint, yaxis = L"\sigma_{\omega}", xaxis = L"\sigma_u", legend = :bottomright, markersize = 5, label = ["Differential" "Algebraic"])
Plots.plot!([0,1], [0,1], line = :dash, label = "", color = colorant"gray0")
png(joinpath(@__DIR__,"../plots/ensemble/short_circ/surv_u_omega"))


@df df groupedhist(:SNBS, group = :Constraint, xlims = [0, 1], bins = 12, legend = :topleft, fillalpha = 0.7, color = [colorant"lightseagreen" colorant"firebrick1" ],bar_position = :stack,   xaxis = L"\beta_{\omega, u}", bar_edges = false, linewidth = 0, label = ["Differential" "Algebraic"])
@df df groupedhist(:SURV, group = :Constraint, xlims = [0, 1], bins = 12, legend = :topleft, fillalpha = 0.7, color = [colorant"deeppink" colorant"green" ],bar_position = :stack,  xaxis = L"\sigma_{\omega, u}", bar_edges = false, linewidth = 0, label = ["Differential" "Algebraic"])


#########################################
#           Infeasabilities             #
#########################################

get_kullback_leibler_divergence(df_cons[!, :SNBS], 1 .- df_cons[!, :INFEASIBLE], 0:0.001:1)
get_kullback_leibler_divergence(df_dyn[!, :SNBS], 1 .- df_dyn[!, :INFEASIBLE], 0:0.001:1)
get_kullback_leibler_divergence(df_cons[!, :SURV], 1 .- df_cons[!, :INFEASIBLE], 0:0.001:1)
get_kullback_leibler_divergence(df_dyn[!, :SURV], 1 .- df_dyn[!, :INFEASIBLE], 0:0.001:1)

jensen_shannon_divergence(df_cons[!, :SNBS], 1 .- df_cons[!, :INFEASIBLE], 0:0.001:1)
jensen_shannon_divergence(df_dyn[!, :SNBS], 1 .- df_dyn[!, :INFEASIBLE], 0:0.001:1)
jensen_shannon_divergence(df_cons[!, :SURV], 1 .- df_cons[!, :INFEASIBLE], 0:0.001:1)
jensen_shannon_divergence(df_dyn[!, :SURV], 1 .- df_dyn[!, :INFEASIBLE], 0:0.001:1)




using KernelDensity
using StatsPlots

p_xy = kde(df_dyn[!, :SNBS], 1 .- df_dyn[!, :INFEASIBLE])

plot(p_xy, xlims = [0.85, 1], ylims = [0.85, 1])

function mutual_information(data_1, data_2)
    I_XY = 0.0

    P_XY = kde((data_1, data_2)) # joint probability distribution
    P_X = kde(data_1, P_XY.x) # marginal probability distribution
    P_Y = kde(data_2, P_XY.y) # marginal probability distribution

    X = P_X.x
    Y = P_Y.x
    
    P_X = P_X.density
    P_Y = P_Y.density
    P_XY = P_XY.density

    for y in 1:length(Y)
        for x in 1:length(X)
            I_XY += P_XY[x, y] * log(P_XY[x, y] / (P_X[x] * P_Y[y]))
            if isnan(I_XY)
                return x, y, P_XY, P_X, P_Y
            end
        end
    end
    return I_XY
end



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