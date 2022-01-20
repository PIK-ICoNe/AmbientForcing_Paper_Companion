module BS_DAE
    using CSV
    using DataFrames
    using PowerDynamics
    using OrdinaryDiffEq
    using ForwardDiff
    using LinearAlgebra
    using SlurmClusterManager
    using Statistics
    using SyntheticNetworks
    using TreeNodeClassification
    using AmbientForcing
    using Distributions
    using LaTeXStrings
    using Plots.Measures
    using GraphMakie
    using Plots, StatsPlots, Statistics
    using LightGraphs, GraphMakie
    using Colors, ColorSchemes
    using Distributed
    using DelimitedFiles
    using OrderedCollections
    using KernelDensity
    using StatsBase
    using Dierckx

    # load power grid library
    include(joinpath("powergrids/", "PhaseAmplitudeOscillator.jl"))
    include(joinpath("powergrids/", "ieee-rts-96-master/PowerDynamics_dyn.jl"))
    export SchifferApprox, get_ieee_96

    # load utils
    include("utils.jl")
    export find_power_grid, half_success_half_failure, perturbed_nodes, node_to_var_idx,  constraint_vector, group_consecutive
    export filter_degree_Ax, filter_degree, sort_by_distance, set_perturbation, check_eigenvalues, snbs_fault_statistics
    export kullback_leibler_divergence, jensen_shannon_divergence
    
    # load plotting functions
    include("plotting.jl")
    export plot_vs_only_dyn, plot_vs, plot_res, spread_violin_plot, my_graph_plot
    export density_degree_plot, plot_fault_simulation, density_class_plot, density_class_plot!

    # load code SNBS simulations
    include("SNBS.jl")
    export simulate_perturbation, pd_node_idx, random_perturbation, parallel_snbs_surv

    include("spreadability.jl")
    export spread_parallel

    include("pg_ensemble.jl")
    export ensemble_calc
end # module
