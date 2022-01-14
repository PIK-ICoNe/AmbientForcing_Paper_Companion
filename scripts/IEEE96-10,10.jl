using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
Pkg.instantiate()
using DataFrames, CSV
using TreeNodeClassification, DelimitedFiles
using BS_DAE

## IEEE 96 Test Case


pg_ieee96 = get_ieee_96()
cons_idx, dyn_idx, constraint_vec = constraint_vector(pg_ieee96)

# Load Simulated Data
ieee_96_all = readdlm(joinpath(@__DIR__,"../data/SNBS_SURV_All_Nodes_73_SampleSize_500_TolBS_0.1_SimTime_1000.0_TauMax_1.0_Dist_[[-0.5, 0.5], [-10.0, 10.0]].txt")) 
num_tries = 500

# Use half failure half success to estimate the standard error
β_96, eβ_96 = half_success_half_failure(ieee_96_all[:,1], num_tries)
β_96_u, eβ_96_u = half_success_half_failure(ieee_96_all[:,2], num_tries)
β_96_ω, eβ_96_ω = half_success_half_failure(ieee_96_all[:,3], num_tries)
σ_96, eσ_96 = half_success_half_failure(ieee_96_all[:,4], num_tries)
σ_96_u, eσ_96_u = half_success_half_failure(ieee_96_all[:,5], num_tries)
σ_96_ω, eσ_96_ω = half_success_half_failure(ieee_96_all[:,6], num_tries)

plot_vs_only_dyn(pg_ieee96, β_96_u, β_96, eβ_96_u, eβ_96, L"\beta_{u}", L"\beta", colorant"mediumorchid", [0.5, 1])
png(joinpath(@__DIR__,"../plots/ieee96/-10,10/snbs_96_all_vs_u"))

plot_vs_only_dyn(pg_ieee96, β_96_u, β_96_ω , eβ_96_u, eβ_96, L"\beta_{u}", L"\beta_{\omega}", colorant"coral1", [0.5, 1])
png(joinpath(@__DIR__,"../plots/ieee96/-10,10/snbs_96_omgea_vs_u"))

plot_vs_only_dyn(pg_ieee96, σ_96_u, σ_96, eσ_96_u, eσ_96, L"\sigma_u", L"\sigma", colorant"steelblue", [0.50, 1])
png(joinpath(@__DIR__,"../plots/ieee96/-10,10/surv_96_all_vs_u"))

plot_vs_only_dyn(pg_ieee96, σ_96_u, σ_96_ω, eσ_96_u, eσ_96, L"\sigma_u", L"\sigma_{\omega}",  colorant"teal", [0.5, 1])
png(joinpath(@__DIR__,"../plots/ieee96/-10,10/surv_96_omgea_vs_u"))

df = DataFrame(Constraint = constraint_vec, SNBS = β_96, SURV = σ_96)

@df df groupedhist(:SNBS, group = :Constraint, xlims = [0.4, 1], fillalpha = 0.8, bins = 12,legend = :topleft, color = [colorant"coral1" colorant"firebrick1" ],bar_position = :stack,  xaxis = L"\beta_{\omega, u}", bar_edges = false, linewidth = 0, label = ["Differential" "Algebraic"])
png(joinpath(@__DIR__,"../plots/ieee96/-10,10/snbs_96_histogram"))

@df df groupedhist(:SURV, group = :Constraint, xlims = [0.4, 1], fillalpha= 0.8, bins = 12,legend = :topleft, color = [colorant"teal" colorant"skyblue2"],bar_position = :stack,  xaxis = L"\sigma_{\omega, u}", bar_edges = false, linewidth = 0, label = ["Differential" "Algebraic"])
png(joinpath(@__DIR__,"../plots/ieee96/-10,10/surv_96_histogram"))

plot_vs(pg_ieee96, β_96, σ_96, eβ_96, eσ_96, L"\beta_{u, \omega}", L"\sigma_{u,\omega}", colorant"orangered",colorant"darkorchid3", :topleft, [0.4, 1])
png(joinpath(@__DIR__,"../plots/ieee96/-10,10/snbs_surv_96_all"))

####
# Some Graph plots
using CairoMakie
using GraphMakie
using Graphs, FileIO

g = pg_ieee96.graph
path = joinpath(@__DIR__, "../plots/ieee96/-10,10/")

node_color_snbs = fill(colorant"firebrick1", 73)
node_color_snbs[dyn_idx] .= colorant"coral1"

node_color_surv = fill(colorant"teal", 73)
node_color_surv[dyn_idx] .= colorant"steelblue"

g = pg_ieee96.graph
f, ax, p = graphplot(g, node_size = (β_96 .- minimum(σ_96)) * 30, node_color = node_color_snbs)
hidedecorations!(ax); hidespines!(ax)
ax.aspect = DataAspect()
FileIO.save(path * "graph_snbs.png", f)


f, ax, p = graphplot(g, node_size = (σ_96 .- minimum(σ_96)) * 40, node_color = node_color_surv)
hidedecorations!(ax); hidespines!(ax)
ax.aspect = DataAspect()
FileIO.save(path * "graph_surv.png", f)