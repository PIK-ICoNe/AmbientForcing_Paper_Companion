using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
Pkg.instantiate()
using DataFrames, CSV
using TreeNodeClassification, DelimitedFiles
using BS_DAE
default(grid = false, foreground_color_legend = nothing, bar_edges = false, framestyle =:box, msc = :auto, dpi=300, legendfontsize = 11, labelfontsize = 12, tickfontsize = 10)

## IEEE 96 Test Case
pg_ieee96 = get_ieee_96()

cons_idx, dyn_idx, constraint_vec = constraint_vector(pg_ieee96)

# Load Simulated Data
ieee_96_all = readdlm(joinpath(@__DIR__,"../data/SNBS_SURV_All_Nodes_73_SampleSize_500_TolBS_0.1_SimTime_1000.0_TauMax_1.0_Dist_[[-0.1, 0.1], [-3.0, 2.0]].txt")) 
num_tries = 500

# Use half failure half success to estimate the standard error
β_96, eβ_96 = half_success_half_failure(ieee_96_all[:,1], num_tries)
β_96_u, eβ_96_u = half_success_half_failure(ieee_96_all[:,2], num_tries)
β_96_ω, eβ_96_ω = half_success_half_failure(ieee_96_all[:,3], num_tries)
σ_96, eσ_96 = half_success_half_failure(ieee_96_all[:,4], num_tries)
σ_96_u, eσ_96_u = half_success_half_failure(ieee_96_all[:,5], num_tries)
σ_96_ω, eσ_96_ω = half_success_half_failure(ieee_96_all[:,6], num_tries)

findall(β_96 .< 0.999)
findall(σ_96 .< 0.999)

findmax(eβ_96)
findmax(eσ_96)
eβ_96[70] * 100
eσ_96[16] * 100

mean(β_96)
mean(σ_96)

plot_vs_only_dyn(pg_ieee96, β_96_u, β_96, eβ_96_u, eβ_96, L"\beta_{u}", L"\beta", colorant"mediumorchid", [0.5, 1])
plot_vs_only_dyn(pg_ieee96, β_96_u, β_96_ω , eβ_96_u, eβ_96, L"\beta_{u}", L"\beta_{\omega}", colorant"coral1", [0.5, 1])
plot_vs_only_dyn(pg_ieee96, σ_96_u, σ_96, eσ_96_u, eσ_96, L"\sigma_u", L"\sigma", colorant"steelblue", [0.50, 1])
plot_vs_only_dyn(pg_ieee96, σ_96_u, σ_96_ω, eσ_96_u, eσ_96, L"\sigma_u", L"\sigma_{\omega}",  colorant"teal", [0.5, 1])

df = DataFrame(Constraint = constraint_vec, SNBS = β_96, SURV = σ_96)

@df df groupedhist(:SNBS, group = :Constraint, xlims = [0.94, 1], fillalpha = 0.8, bins = 12,legend = :topleft, color = [colorant"coral1" colorant"firebrick1" ],bar_position = :stack,  xaxis = L"\beta_{\omega, u}", bar_edges = false, linewidth = 0, label = ["Differential" "Algebraic"])

@df df groupedhist(:SURV, group = :Constraint, xlims = [0.94, 1], fillalpha= 0.8, bins = 12,legend = :topleft, color = [colorant"teal" colorant"skyblue2"],bar_position = :stack,  xaxis = L"\sigma_{\omega, u}", bar_edges = false, linewidth = 0, label = ["Differential" "Algebraic"])

plot_vs(pg_ieee96, β_96, σ_96, eβ_96, eσ_96, L"\beta_{u, \omega}", L"\sigma_{u,\omega}", colorant"orangered",colorant"darkorchid3", :topleft, [0.94, 1])
