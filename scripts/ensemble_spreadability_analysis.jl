using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
using BS_DAE
using CSV, DataFrames
using KernelDensity, StatsBase, Dierckx
default(grid = false, foreground_color_legend = nothing, bar_edges = false, framestyle =:box, msc = :auto, dpi=300, legendfontsize = 11, labelfontsize = 12, tickfontsize = 10)

input = joinpath(@__DIR__,"../data/Pg_ensemble_method_All_Networks_100_Nodes_100_SampleSize_250_TauMax1.0_Dist[[-2.0, 2.0], [-5.0, 5.0]].csv")
df = DataFrame(CSV.File(input))

df_cons = filter("Constraint" => x -> x == 1, df)
df_dyn = filter("Constraint" => x -> x == 0, df)

#############################
#    Spreadability Plots    #
#############################
# Figure 1 in the Paper
nan_idx = findall(map(x -> isnan(x), df[!, :Abs_Spread]))
delete!(df, nan_idx)

spread_by_dist, spread_sorted, distances_sorted = sort_by_distance(df[!, "Abs_Spread"], df[!, "Constraint_Distance"])
dist_zero =  df[!, "Constraint_Distance"] .* abs.(1 .- df[!, "Constraint"])
spread_by_dist_0, spread_sorted_0, distances_sorted_0 = sort_by_distance(df[!, "Abs_Spread"], dist_zero)

plt1 = spread_violin_plot(distances_sorted, spread_sorted, colorant"orange2", colorant"steelblue", latexstring(string(raw"\min (d", "_{cn})")))
plt2 = spread_violin_plot(distances_sorted_0, spread_sorted_0, colorant"peachpuff3", colorant"royalblue1", latexstring(string(raw"\min (d", "_{c})")))

###############################
#  Amplification vs. Degree   #
###############################
# Figure in the Appendix

Ax_degree, Ax_d_δ = filter_degree_Ax(df)
my_range = 0:0.01:4

dens = kde(df[!,:AX])
dens = Spline1D(dens.x, dens.density)

Plots.plot(my_range, dens(my_range) ./ sum(dens(my_range)),  lw = 3, xaxis = L"A", yaxis = L"p(A)", legend = false)
Plots.plot(string.(collect(1:length(Ax_degree))), Ax_degree, ribbon = (Ax_d_δ[:,1], Ax_d_δ[:,2]), fillalpha = 0.3, lw = 3, xaxis = L"d", yaxis = L"A", legend = false)