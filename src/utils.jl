"""
    sort_by_distance(data, distance)

Sorts an array of data by a distance measure.
"""
function sort_by_distance(data, distance)
    distances_sorted, data_sorted = getindex.((distance, data), (sortperm(distance),))
    max_dist = Int(distances_sorted[end])
    data_by_dist = Array{Array}(undef, max_dist)

    for i in 1:max_dist
        idx = findall(x -> x == i, distances_sorted)
        data_by_dist[i] = data_sorted[idx]
    end
    return data_by_dist, data_sorted, distances_sorted
end

"""
    group_consecutive(b)
"""
function group_consecutive(b)
    x = b[2:end] .- b[1:end-1] .== 1 

    chunks = []
    temp = 1
    for i in 1:length(x)
        if x[i] == 0
            append!(chunks, [b[temp:i]])
            temp = i + 1
        end
        if i == length(x)
            append!(chunks, [b[temp:i+1]])
        end
    end
    return chunks
end

"""
    filter_degree(df::DataFrame)

Calculates the mean Snbs and Surv for each degree.
"""
function filter_degree(df::DataFrame)
    max_degree = maximum(df[!, "degree"])

    nodes_per_degree = map(i -> nrow(filter("degree" => x -> x == i, df)), 1:max_degree)
    last_degree = findlast(nodes_per_degree .> 10)

    SNBS_degree = zeros(last_degree)
    SURV_degree = zeros(last_degree)
    SURV_δ = zeros(last_degree, 2)
    SNBS_δ = zeros(last_degree, 2)

    for i in 1:last_degree
        my_df = filter("degree" => x -> x == i, df)

        SURV_degree[i] = mean(my_df[!, "SURV"])
        SNBS_degree[i] = mean(my_df[!, "SNBS"])

        SURV_δ[i, 1] =  abs(percentile(my_df[!, "SURV"], 15.9) - SURV_degree[i]) 
        SURV_δ[i, 2] =  abs(percentile(my_df[!, "SURV"], 84.1) - SURV_degree[i])

        SNBS_δ[i, 1] =  abs(percentile(my_df[!, "SNBS"], 15.9) - SNBS_degree[i])
        SNBS_δ[i, 2] =  abs(percentile(my_df[!, "SNBS"], 84.1) - SNBS_degree[i])
    end
    return SNBS_degree, SNBS_δ, SURV_degree, SURV_δ
end

"""
    filter_degree_Ax(df::DataFrame)

Calculates the mean Amplification for each degree.
"""
function filter_degree_Ax(df::DataFrame)
    max_degree = maximum(df[!, "degree"])

    nodes_per_degree = map(i -> nrow(filter("degree" => x -> x == i, df)), 1:max_degree)
    last_degree = findlast(nodes_per_degree .> 10)

    AX_degree = zeros(last_degree)
    AX_δ = zeros(last_degree, 2)

    for i in 1:last_degree
        my_df = filter("degree" => x -> x == i, df)

        AX_degree[i] = mean(my_df[!, "AX"])

        AX_δ[i, 1] =  abs(percentile(my_df[!, "AX"], 15.9) - AX_degree[i]) 
        AX_δ[i, 2] =  abs(percentile(my_df[!, "AX"], 84.1) - AX_degree[i])
    end
    return AX_degree, AX_δ
end


"""
    constraint_vector(pg::PowerGrid)

Returns the indexes of the constraint and dynamical nodes as a vector.
"""
function constraint_vector(pg::PowerGrid)
    rpg = rhs(pg)
    num_nodes = length(pg.nodes)
    M = rpg.mass_matrix
    constraint_vector = zeros(num_nodes)
    for node in 1:num_nodes
        node_vars = node_to_var_idx(rpg, node)
        if M[node_vars[1], node_vars[1]] == 0 # M[i,i] = 0 means constraint
            constraint_vector[node] = 1
        end
    end
    dyn_idx = findall(constraint_vector .== 0.0)
    cons_idx = findall(abs.(constraint_vector .- 1 ) .== 0.0)
    return cons_idx, dyn_idx, constraint_vector 
end

"""
    node_to_var_idx(f::ODEFunction, node::Int64)

Takes a ODEFunction from PowerDynamics and returns the indexes of the dynamical
variables associated to a node.
"""
function node_to_var_idx(f::ODEFunction, node::Int64)
    M = f.mass_matrix
    len_M = size(M, 1)
    var_array = []

    for i in 1:len_M
        str = string(f.syms[i])
        idx = findlast('_', str)
        node_num = parse(Int64, str[idx + 1:end])
        
        if node_num > node
            return var_array
        end

        if node_num == node
            append!(var_array, i)
        end
    end
    return var_array
end

"""
    perturbed_nodes(rhs, op, sol)

Returns the node indices at which the solution deviates from the initial condition.
"""
function perturbed_nodes(rhs, op, sol)
    dev_idx = findall(map(x -> x > 0.01, abs.(op.vec - sol(0.1))))
    str = string.(rhs.syms[dev_idx])
    idx = findlast.('_', str)
    disturbed_nodes = []

    for i in 1:length(idx)
        my_str = str[i]
        node_num = parse(Int64, my_str[idx[i] + 1:end])
        if node_num ∈ disturbed_nodes
            continue
        else
            append!(disturbed_nodes, node_num)
        end
    end
    return disturbed_nodes
end

"""
    half_success_half_failure(x, L::Int)

"""
function half_success_half_failure(x, L::Int)
    s = x * L
    x_tilde = (s .+ 0.5) ./ (L + 1)
    e_tilde = sqrt.(x_tilde .* (1 .- x_tilde) ./ (L + 1))
    return x_tilde, e_tilde
end

"""
    set_perturbation(pg::PowerGrid, rpg::ODEFunction, node::Int, op, pert_vec)

Returns the specified perturbation using Ambient Forcing
"""
function set_perturbation(pg::PowerGrid, rpg::ODEFunction, node::Int, op, pert_vec)
    idx = pd_node_idx(pg, rpg, node, "All")
    h = zeros(length(op.vec))

    if length(idx) == 2
        h[idx] .= pert_vec[1:2]
    elseif length(idx) == 3
        h[idx] .= pert_vec[1:3]
    elseif length(idx) == 4
        h[idx] .= pert_vec[1:4]
    end
    ambient_forcing(rpg, op.vec, 1.0, h) # solve the differential equation
end

"""
    check_eigenvalues(pg::PowerGrid, s::State)

Calculates the eigenvalues λ of the jacobian of the right hand side of
the power grid at the state s.
The sign of the real part of the eigenvalues will decide whether a fixed point
will be attracting (λ_max < 0) or repelling (λ_max > 0).
"""
function check_eigenvalues(pg::PowerGrid, s::State)
    rpg = rhs(pg)
    M = Array(rpg.mass_matrix) # ohne slack failt
    f!(dx, x) = rpg(dx, x, nothing, 0.)
    j(x) = (dx = similar(x); ForwardDiff.jacobian(f!, dx, x))
    λ = eigvals(j(s.vec) * pinv(M) * M) .|> real |> extrema
    stable = isapprox(last(λ), 0, atol=1e-8)
    return λ, stable
end

"""
    snbs_fault_statistics(df::DataFrame)

Calculate the percentage of nodes where a certain fault occurs.
"""
function snbs_fault_statistics(df::DataFrame)
    # all nodes where a voltage appears
    vol_drop_nodes = findall((1 .- df[!, :INFEASIBLE]) .- df[!, :SNBS_u] .> 0.001) 
    p_vol = (length(vol_drop_nodes) / length(df[!, :SNBS])) * 100
    println("Voltage Drops occur at: ", string(p_vol) , " % of nodes")

    # all nodes where a desynronization appears
    desyn_nodes = findall((1 .- df[!, :INFEASIBLE]) .- df[!, :SNBS_ω] .> 0) 
    p_desyn = (length(desyn_nodes) / length(df[!, :SNBS])) * 100
    println("Desynronizations occur at: ", string(p_desyn) , " % of nodes") 

    # all nodes where voltage drops appear and no desynronization
    vol_drop_only = findall(df[!, :SNBS_ω] .- df[!, :SNBS_u] .> 0.001) 
    p_vol_only = (length(vol_drop_only) / length(df[!, :SNBS])) * 100
    println("Voltage Drops without Desynronization occur at: ", string(p_vol_only) , " % of nodes")
   
    # all nodes where infeasible power flows occur
    infeas_nodes = findall(df[!, :INFEASIBLE] .> 0.0) 
    p_infeas = (length(infeas_nodes) / length(df[!, :SNBS])) * 100
    println("Infeasible Solutions occur at: ", string(p_infeas) , " % of nodes")

    return p_vol, p_desyn, p_vol_only, p_infeas
end

"""
    find_power_grid(pg_idx::Int)

Load the power grid from the data
"""
function find_power_grid(pg_idx::Int)
    # load the power grid and all related stuff
    filename = joinpath(@__DIR__,"../data/Ensemble/pg_method_Voltage_pg_idx_" * string(pg_idx))
    pg = read_powergrid(filename, Json)
    return pg, rhs(pg), find_operationpoint(pg)
end

"""
    kullback_leibler_divergence(P, Q, X)

The Kullback–Leibler divergence is a statistical distance: a measure of how one probability distribution Q is different from a second, reference probability distribution P.
- `P`: Reference probability distribution
- `Q`: Probability distribution
- `X`: Interval of possible outcomes to run over (in Basin Stability [0,1])

https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence
"""
function kullback_leibler_divergence(P, Q, X)
    D_KL = 0
    for x in X
        D_KL += abs(P(x)) * log(abs(P(x))/ abs(Q(x)))
    end
    return D_KL
end

"""
    jensen_shannon_divergence(data_1, data_2, X)

Calculates the Jensen-Shannon-divergence of the probability densities given by data_1 and data_2.
The Jensen-Shannon divergence, is a symmetric and bounded version of the Kullback-Leibler divergence.
- `data_1`: Reference data set
- `data_2`: Data set
- `X`: Interval of possible outcomes to run over (in Basin Stability [0,1])
"""
function jensen_shannon_divergence(data_1, data_2, X)
    P = kde(data_1)
    Q = kde(data_2)

    P = Spline1D(P.x, P.density ./ sum(P.density)) 
    Q = Spline1D(Q.x, Q.density ./ sum(Q.density))
    M = 1/2 .* (P(X) + Q(X))
    M = Spline1D(X, M)

    D_PM = kullback_leibler_divergence(P, M, X)
    D_QM = kullback_leibler_divergence(Q, M, X)

    JSD = 1/2 * (D_PM + D_QM)

    return JSD
end

import PowerDynamics: convert_node
using PowerDynamics:_map_complex
export convert_node
function convert_node(node)
    name = get(node, "name", nothing)
    type = get(node, "type", nothing)
    params = get(node, "params", [])
    sym_params = Dict(Symbol(k) => _map_complex(v) for (k, v) in params)
    if type == "SwingEq"
        node = SwingEq(;sym_params...)
    elseif type == "SwingEqLVS"
        node = SwingEqLVS(;sym_params...)
    elseif type == "FourthOrderEq"
        node = FourthOrderEq(;sym_params...)
    elseif type == "FourthOrderEqGovernorExciterAVR"
        node = FourthOrderEqGovernorExciterAVR(;sym_params...)
    elseif type == "SlackAlgebraic"
        node = SlackAlgebraic(;sym_params...)
    elseif type == "PQAlgebraic"
        node = PQAlgebraic(;sym_params...)
    elseif type == "PVAlgebraic"
        node = PVAlgebraic(;sym_params...)
    elseif type == "PVAlgebraic"
        node = PVAlgebraic(;sym_params...)
    elseif type == "VSIMinimal"
        node = VSIMinimal(;sym_params...)
    elseif type == "VSIVoltagePT1"
        node = VSIVoltagePT1(;sym_params...)
    elseif type == "CSIMinimal"
        node = CSIMinimal(;sym_params...)
    elseif type == "ExponentialRecoveryLoad"
        node = ExponentialRecoveryLoad(;sym_params...)
    elseif type == "VoltageDependentLoad"
        node = VoltageDependentLoad(;sym_params...)
    elseif type == "SchifferApprox"
        node = SchifferApprox(;sym_params...)
    else
        throw(ArgumentError("Invalid type: $type"))
    end
    if typeof(name) == Nothing
        node
    else
        (name,node)
    end
end