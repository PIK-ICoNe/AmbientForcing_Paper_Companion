using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
Pkg.instantiate()
using DataFrames, CSV, DelimitedFiles
using KernelDensity
using BS_DAE

num_tries = 500
cons_idx, dyn_idx, constraint_vec = constraint_vector(pg_ieee96)

# Load simulated Data
#vol_01 = readdlm(joinpath(@__DIR__,"../data/Voltage Sweep/SNBS_SURV_Voltage_Nodes_73_SampleSize_500_TolBS_0.1_SimTime_750.0_TauMax_1.0_Dist_[[-0.1, 0.1], [-3.0, 2.0]].txt")) 
#vol_05 = readdlm(joinpath(@__DIR__,"../data/Voltage Sweep/SNBS_SURV_Voltage_Nodes_73_SampleSize_500_TolBS_0.1_SimTime_750.0_TauMax_1.0_Dist_[[-0.5, 0.5], [-3.0, 2.0]].txt")) 
vol_075 = readdlm(joinpath(@__DIR__,"../data/Voltage Sweep/SNBS_SURV_Voltage_Nodes_73_SampleSize_500_TolBS_0.1_SimTime_750.0_TauMax_1.0_Dist_[[-0.75, 0.75], [-3.0, 2.0]].txt")) 
vol_1 = readdlm(joinpath(@__DIR__,"../data/Voltage Sweep/SNBS_SURV_Voltage_Nodes_73_SampleSize_500_TolBS_0.1_SimTime_750.0_TauMax_1.0_Dist_[[-1.0, 1.0], [-3.0, 2.0]].txt")) 
vol_15 = readdlm(joinpath(@__DIR__,"../data/Voltage Sweep/SNBS_SURV_Voltage_Nodes_73_SampleSize_500_TolBS_0.1_SimTime_750.0_TauMax_1.0_Dist_[[-1.5, 1.5], [-3.0, 2.0]].txt")) 
vol_2 = readdlm(joinpath(@__DIR__,"../data/Voltage Sweep/SNBS_SURV_Voltage_Nodes_73_SampleSize_500_TolBS_0.1_SimTime_750.0_TauMax_1.0_Dist_[[-2.0, 2.0], [-3.0, 2.0]].txt")) 
vol_3 = readdlm(joinpath(@__DIR__,"../data/Voltage Sweep/SNBS_SURV_Voltage_Nodes_73_SampleSize_500_TolBS_0.1_SimTime_750.0_TauMax_1.0_Dist_[[-3.0, 3.0], [-3.0, 2.0]].txt")) 
vol_4 = readdlm(joinpath(@__DIR__,"../data/Voltage Sweep/SNBS_SURV_Voltage_Nodes_73_SampleSize_500_TolBS_0.1_SimTime_750.0_TauMax_1.0_Dist_[[-4.0, 4.0], [-3.0, 2.0]].txt")) 
vol_5 = readdlm(joinpath(@__DIR__,"../data/Voltage Sweep/SNBS_SURV_Voltage_Nodes_73_SampleSize_500_TolBS_0.1_SimTime_750.0_TauMax_1.0_Dist_[[-5.0, 5.0], [-3.0, 2.0]].txt")) 


#β_01, eβ_01 = half_success_half_failure(vol_01[:,1], num_tries)
#β_05, eβ_05 = half_success_half_failure(vol_05[:,1], num_tries)
β_075, eβ_075 = half_success_half_failure(vol_075[:,1], num_tries)
β_1, eβ_1 = half_success_half_failure(vol_1[:,1], num_tries)
β_15, eβ_15 = half_success_half_failure(vol_15[:,1], num_tries)
β_2, eβ_2 = half_success_half_failure(vol_2[:,1], num_tries)
β_3, eβ_3 = half_success_half_failure(vol_3[:,1], num_tries)
β_4, eβ_4 = half_success_half_failure(vol_4[:,1], num_tries)
β_5, eβ_5 = half_success_half_failure(vol_5[:,1], num_tries)
β_5, eβ_5 = half_success_half_failure(vol_5[:,1], num_tries)

#σ_01, eσ_01 = half_success_half_failure(vol_01[:,4], num_tries)
#σ_05, eσ_05 = half_success_half_failure(vol_05[:,4], num_tries)
σ_075, eσ_075 = half_success_half_failure(vol_075[:,4], num_tries)
σ_1, eσ_1 = half_success_half_failure(vol_1[:,4], num_tries)
σ_15, eσ_15 = half_success_half_failure(vol_15[:,4], num_tries)
σ_2, eσ_2 = half_success_half_failure(vol_2[:,4], num_tries)
σ_3, eσ_3 = half_success_half_failure(vol_3[:,4], num_tries)
σ_4, eσ_4 = half_success_half_failure(vol_4[:,4], num_tries)
σ_5, eσ_5 = half_success_half_failure(vol_5[:,4], num_tries)

f_075, ef_075 = half_success_half_failure(vol_075[:,6], num_tries)
f_1, ef_1 = half_success_half_failure(vol_1[:,6], num_tries)
f_15, ef_15 = half_success_half_failure(vol_15[:,1], num_tries)
f_2, ef_2 = half_success_half_failure(vol_2[:,6], num_tries)
f_3, ef_3 = half_success_half_failure(vol_3[:,6], num_tries)
f_4, ef_4 = half_success_half_failure(vol_4[:,6], num_tries)
f_5, ef_5 = half_success_half_failure(vol_5[:,6], num_tries)
f_5, ef_5 = half_success_half_failure(vol_5[:,6], num_tries)

lab = [L"\pm 0.75" L"\pm 1.0" L"\pm 1.5" L"\pm 2.0" L"\pm 3.0" L"\pm 4.0" L"\pm 5.0"]

Plots.violin(lab, [β_075[dyn_idx] β_1[dyn_idx] β_15[dyn_idx] β_2[dyn_idx] β_3[dyn_idx] β_4[dyn_idx] β_5[dyn_idx]], side = :left, color = colorant"coral1", linewidth=0, label=["Differential" "" "" "" "" "" ""])
Plots.violin!(lab, [β_075[cons_idx] β_1[cons_idx] β_15[cons_idx] β_2[cons_idx] β_3[cons_idx] β_4[cons_idx] β_5[cons_idx]], ylabel = L"\beta" ,color = colorant"firebrick1",side=:right, linewidth=0, label=["Algebraic" "" "" "" "" "" ""], legend = :bottomleft)
png(joinpath(@__DIR__,"../plots/ieee96/snbs_voltage_sweep_violin"))


Plots.violin(lab, [σ_075[cons_idx] σ_1[cons_idx] σ_15[cons_idx] σ_2[cons_idx] σ_3[cons_idx] σ_4[cons_idx] σ_5[cons_idx]], ylabel = L"\sigma", color = colorant"steelblue",side = :right, linewidth = 0, label=["Algebraic" "" "" "" "" "" ""], legend = :bottomleft)
Plots.violin!(lab, [σ_075[dyn_idx] σ_1[dyn_idx] σ_15[dyn_idx] σ_2[dyn_idx] σ_3[dyn_idx] σ_4[dyn_idx] σ_5[dyn_idx]], side=:left, color = colorant"teal", linewidth = 0, label=["Differential" "" "" "" "" "" ""])
png(joinpath(@__DIR__,"../plots/ieee96/surv_voltage_sweep_violin"))

#########################################
#           Infeasabilities             #
#########################################

1 .- f_5

abs.(β_5 .- (1 .- f_5))

maximum(β_15 .- f_15)

StatsPlots.density(β_075 .- (1 .-f_075), lw = 3, label = "Algebraic")
StatsPlots.density!(β_1 .- (1 .- f_1), lw = 3, label = "Algebraic")
StatsPlots.density!(abs.(β_15 .- (f_15)), lw = 3, label = "Algebraic")
StatsPlots.density(β_2 .- (1 .- f_2), lw = 3, label = "Algebraic")
StatsPlots.density!(β_3 .- (1 .- f_3), lw = 3, label = "Algebraic")
StatsPlots.density!(β_4 .- (1 .- f_4), lw = 3, label = "Algebraic")
StatsPlots.density!(β_5 .- (1 .- f_5), lw = 3, label = "Algebraic", legend = false)

scatter(1 .- f_2, β_2)
Plots.plot!([0,1], [0,1], line = :dash, label = "", color = colorant"gray0")
