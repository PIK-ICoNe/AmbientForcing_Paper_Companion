"""
    random_PD_grid(N::Int)

Generates a random power grid using SyntheticNetworks and then turns it into a PowerDynamics.PowerGrid type.
"""
function random_PD_grid(N::Int; τ_P = 5.0, τ_Q = 8.0, K_P = 5, K_Q = 0.1)
    while true
        # Generates the random power grid graph using SyntheticNetworks
        RPG = RandomPowerGrid(N, 1, 1/5, 3/10, 1/3, 1/10, 0.0)
        pg = generate_graph(RPG)
        e = edges(pg.graph)

        # collects all edges and turns them into StaticLines
        from_vec = src.(e)
        to_vec = dst.(e)

        lines = Array{Any}(undef, length(e))
        nodes = Array{Any}(undef, nv(pg.graph))

        P_vec = rand(Uniform(-1, 1), N)

        for l in 1:length(e)
            lines[l] = StaticLine(from = from_vec[l], to = to_vec[l], Y = -3im) 
        end
        nodes[1] = SlackAlgebraic(U = complex(1.0))

        for n in 2:nv(pg.graph)
            nodes[n] = PVAlgebraic(P = P_vec[n], V = 1.0)
        end
        pg_cons = PowerGrid(nodes, lines)
        op_cons = find_operationpoint(pg_cons, sol_method=:rootfind)#, solve_powerflow = true)

        nodes = Array{Any}(undef, nv(pg.graph))
        nodes[1] = SlackAlgebraic(U = complex(1.0))

        # Randomly turns nodes into PQAlgebraic or PhaseAmplitudeOscillators
        for n in 2:nv(pg.graph)
            α = rand()
            if α > 0.5 # Grid Forming inverter 
                nodes[n] = SchifferApprox(τ_P = τ_P, τ_Q = τ_Q, K_P = K_P, K_Q = K_Q, V_r = 1.0, P = op_cons[n, :p], Q = op_cons[n, :q], Y_n = 0)
            else # Grid Following Inverters
                nodes[n] = PQAlgebraic(P = op_cons[n, :p], Q = op_cons[n, :q]) 
            end
        end

        pg = PowerGrid(nodes, lines)
        op = find_operationpoint(pg, sol_method=:rootfind)#, solve_powerflow = true)

        if all(isapprox.(op[:, :v], 1.0))
            _, stable = check_eigenvalues(pg, op) # check if operation point is linearly stable
            if stable 
                sol = simulate_perturbation(rhs(pg), op.vec, (0.0, 25.0))
                if sol.retcode == :Success
                    return pg, op
                end
            end
        end
    end
end 

"""
    spread_ensemble(dist_args, tau_max::Float64, num_networks::Int64, throws::Int64, num_nodes::Int64)

Calculates the relevant measures for an ensemble of randomly created synthetical power grids. 
"""
function ensemble_calc(dist_vec, tau_max::Float64, num_networks::Int64, sample_size::Int64, num_nodes::Int64, method::String; τ_P = 5.0, τ_Q = 8.0, K_P = 5, K_Q = 0.1, spread_run = false, tol = false, simtime = 500.0)
    if tol != false
        file_name = "Pg_ensemble_method_" * method * "_Networks_" * string(num_networks) * "_Nodes_" * string(num_nodes) * "_SampleSize_" * string(sample_size) * "_TauMax" * string(tau_max) * "_Dist" * string(dist_vec) * "_tol_" * string(tol) * ".csv"
    elseif simtime != 500.0
        file_name = "Pg_ensemble_method_" * method * "_Networks_" * string(num_networks) * "_Nodes_" * string(num_nodes) * "_SampleSize_" * string(sample_size) * "_TauMax" * string(tau_max) * "_Dist" * string(dist_vec) * "_simtime_" * string(simtime) * ".csv"
    elseif τ_P != 5.0 || τ_Q != 8.0 || K_P != 5 || K_Q != 0.1
        file_name = "Pg_ensemble_method_" * method * "_Networks_" * string(num_networks) * "_Nodes_" * string(num_nodes) * "_SampleSize_" * string(sample_size) * "_Dist" * string(dist_vec) * "_" * string(τ_P) *  "_" * string(τ_Q) * "_" * string(K_P) * "_" * string(K_Q) * ".csv"
    else
        file_name = "Pg_ensemble_method_" * method * "_Networks_" * string(num_networks) * "_Nodes_" * string(num_nodes) * "_SampleSize_" * string(sample_size) * "_TauMax" * string(tau_max) * "_Dist" * string(dist_vec) * ".csv"
    end

    output = joinpath(@__DIR__, "../data/", file_name)
    
    for i in 1:num_networks
        # sampling a random powergrid
        pg, op = random_PD_grid(num_nodes, τ_P = τ_P, τ_Q = τ_Q, K_P = K_P, K_Q = K_Q)
        pg_idx = rand(Int)

        power = map(x -> pg.nodes[x].P, 2:num_nodes)
        power = vcat("Slack", power)

        dist_min, cons_vec = min_constraint_distance(pg) # if cons_vec = 1 the node is a constraint
        class_vec = full_node_classification(pg.graph, 1000, 5) # Topological Classes
        
        snbs, snbs_u, snbs_ω, surv, surv_u, surv_ω, infeasible = parallel_snbs_surv(pg, op, simtime, sample_size, dist_vec, 0.1, tau_max, method, false, tol = tol) # Basin Stability and Survivability
        
        if spread_run
            abs_spread, Ax = spread_parallel(pg, op, dist_vec[1:2], tau_max, sample_size) # Spreadability
            df = DataFrame(PG_IDX = fill(pg_idx, num_nodes), Abs_Spread = abs_spread, AX = Ax, Power = power, INFEASIBLE = infeasible, Constraint_Distance = dist_min, Constraint = cons_vec, SNBS = snbs, SNBS_u = snbs_u, SNBS_ω = snbs_ω, SURV = surv, SURV_u = surv_u, SURV_ω = surv_ω, NodeClass = class_vec, degree = degree(pg.graph), betweenness = betweenness_centrality(pg.graph))
        else
            df = DataFrame(PG_IDX = fill(pg_idx, num_nodes), Power = power, INFEASIBLE = infeasible, Constraint_Distance = dist_min, Constraint = cons_vec, SNBS = snbs, SNBS_u = snbs_u, SNBS_ω = snbs_ω, SURV = surv, SURV_u = surv_u, SURV_ω = surv_ω, NodeClass = class_vec, degree = degree(pg.graph), betweenness = betweenness_centrality(pg.graph))
        end
        
        write_powergrid(pg, joinpath(@__DIR__, "../data/Ensemble/pg_method_" * method  * "_pg_idx_" * string(pg_idx)), Json)

        if isfile(output)
            CSV.write(output, df, append = true)
        else
            CSV.write(output, df, append = false)
        end
    end
    return true
end

"""
    min_constraint_distance(pg::PowerGrid)

Calculates the minimal distance to the next constraint in a power grid.
"""
function min_constraint_distance(pg::PowerGrid)
    num_nodes = length(pg.nodes)
    min_distance = Array{Float64}(undef, num_nodes)
    constraint_idx, dynamical_vec, constraint_vec = constraint_vector(pg)
    f = rhs(pg)

    # handling the case of no constraints
    if f.mass_matrix == true * I 
        min_distance .= num_nodes

    # handling the case of only constraints
    elseif length(constraint_idx) == num_nodes 
        min_distance = zeros(num_nodes)
    else
        for node in 1:num_nodes
            cons_vec = copy(constraint_idx) # we do not want to overwrite the original vector

            # getting the shortest path lengths
            edge_dist = bellman_ford_shortest_paths(pg.graph, node)
            dists = edge_dist.dists
                
            if node ∈ constraint_idx
                cons_vec = filter(x -> x != node, cons_vec)
                
                if length(cons_vec) == 0 # handling the case of this node being the only constraint
                    min_distance[node] = num_nodes
                else
                    min_distance[node] = minimum(dists[cons_vec])
                end
            else
                min_distance[node] = minimum(dists[cons_vec])
            end
        end
    end
    return min_distance, constraint_vec
end