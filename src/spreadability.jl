"""
    Spreadability(pg::PowerGrid, op, dist_args, tau_max::Float64, throws::Int64)

Calculates the spreadability of the system. s(z) = ∑ d² * mean(Δz)² / (∑d * ∑ mean(Δz))
"""
function spread(pg::PowerGrid, rpg::ODEFunction, indices, op::State, ur_idx, ui_idx, nodes, dist_args, tau_max::Float64, throws::Int64)
    Δx_abs = zeros(length(pg.nodes), length(pg.nodes))
    for n in 1:length(pg.nodes)
            idx = indices[n]
            dist_vec_node = [sort([(op[n, :u_r]) * dist_args[1][1], (op[n, :u_r]) * dist_args[1][2]]), sort([op[n, :u_i] * dist_args[1][1], op[n, :u_i] * dist_args[1][2]]) ]
            Frand = random_force(rpg, dist_vec_node, Uniform, idx)
            z_new = ambient_forcing(rpg, op.vec, tau_max, Frand) # solve the differential equation

            # getting the displacement after perturbing "node"
            Δx = z_new .- op.vec
            Δx_abs[n, :] = abs.(Δx[ur_idx] .+ 1im * Δx[ui_idx])
    end
    return Δx_abs
end

"""
    spread_parallel(pg::PowerGrid, op, dist_args, tau_max::Float64, sample_size::Int64)

Calculates the Spreadability. Uses pmap and spread() for the parallelization.
"""
function spread_parallel(pg::PowerGrid, op, dist_args, tau_max::Float64, sample_size::Int64)
    rpg, ur_idx, ui_idx, nodes, indices = preliminary_spread(pg)
    
    result = pmap(sample -> spread(pg, rpg, indices, op, ur_idx, ui_idx, nodes, dist_args, tau_max, sample), 1:sample_size)
    Δx_abs = mean(result)

    s_abs = zeros(length(pg.nodes))
    A_x = zeros(length(pg.nodes))

    for n in 1:length(pg.nodes)
        # the voltage of the slack bus is fixed and can not be perturbed! 
        if typeof(pg.nodes[nodes[n]]) == SlackAlgebraic
            s_abs[n] = 0.0
        else
            A_x[n] = sum(Δx_abs[n, :])
            # getting the shortest path length to "node" for all vertices in the network
            edge_dist = bellman_ford_shortest_paths(pg.graph, n)
            dists = edge_dist.dists
            C_dash = length(pg.nodes) / sum(dists) # modified closeness centrality
            # calculating the spreadability
            p_i_x = Δx_abs[n, :] ./ A_x[n]
            s_abs[n] = sum(dists .* p_i_x) * C_dash
        end
    end
    return s_abs, A_x
end

"""
    preliminary_spread(pg::PowerGrid)

Calculates some preliminary things which stay the same during all runs.
"""
function preliminary_spread(pg::PowerGrid)
    rpg = rhs(pg)
    ur_idx = findall(map(variable -> occursin("u_r", variable), string.(rpg.syms)))
    ui_idx = findall(map(variable -> occursin("u_i", variable), string.(rpg.syms)))

    indices = map(n -> pd_node_idx(pg, rpg, n, "Voltage"), 1:length(pg.nodes))

    # handling the type of the pg
    if typeof(pg.nodes) == OrderedDict{String, AbstractNode}
        nodes = collect(keys(pg.nodes)) # we have to give the key to find the node in the dict, but the key can not be used for the numbering
    else 
        nodes = 1:length(pg.nodes) # in Arrays it does not matter
    end
    return rpg, ur_idx, ui_idx, nodes, indices
end