"""
    random_perturbation(idx, rpg, op, dist_vec, tau_max)

Return a random perturbation using the AmbientForcing Algorithm
"""
function random_perturbation(idx, rpg, op, dist_vec, tau_max)
    tau_rand = rand(Uniform(0.0, tau_max))
    Frand = random_force(rpg, dist_vec, Uniform, idx)
    z_new = ambient_forcing(rpg, op.vec, tau_rand, Frand) # solve the differential equation
    return z_new
end

"""
    pd_node_idx(pg::PowerGrid, node::Int, method::String)

Returns index vector of variables at node. Can be used in Random Force
"""
function pd_node_idx(pg::PowerGrid, rpg, node_num::Int, method::String)
    if method == "Voltage"
        str_vec = ["u_r_" * string(node_num), "u_i_" * string(node_num)]
    elseif method == "All"
        if typeof(pg.nodes) == OrderedDict{String, Any}
            node = collect(keys(pg.nodes))[node_num] # we have to give the key to find the node in the dict, but the key can not be used for the numbering
        else 
            node = node_num # in Arrays it does not matter
        end

        if length(symbolsof(pg.nodes[node])) == 2
            str_vec = ["u_r_" * string(node_num), "u_i_" * string(node_num)]
        elseif :ω ∈ symbolsof(pg.nodes[node])
            if :θ ∈ symbolsof(pg.nodes[node])
                str_vec = ["u_r_" * string(node_num), "u_i_" * string(node_num), "ω_" * string(node_num), "θ_" * string(node_num)]
            else
                str_vec = ["u_r_" * string(node_num), "u_i_" * string(node_num), "ω_" * string(node_num)]
            end
        else
            error("Please construct the vector manually")
        end
        else
            error("Please use a valid method.")
    end
    return idx_exclusive(rpg, str_vec)
end

"""
    simulate_perturbation(rpg, x0, tspan)

Simulate State x0 of the function rpg for the time span tspan.
"""
function simulate_perturbation(rpg, x0, tspan; tol = false)
    problem = ODEProblem(rpg, x0, tspan)
    if tol == false
        solution = solve(problem, Rodas4(), force_dtmin = true)
        return solution
    else
        solution = solve(problem, Rodas4(), force_dtmin = true, abstol = tol, reltol = tol)
        return solution
    end
end

"""
    adapted_dist_vec(idx, node::Int, op, dist_vec)

Return dist_vec fitted for the node
"""
function adapted_dist_vec(idx, node::Int, op, dist_vec)
    if length(idx) == 2
        dist_vec_node = [sort([op[node, :u_r] * dist_vec[1][1], op[node, :u_r] * dist_vec[1][2]]), sort([op[node, :u_i] * dist_vec[1][1], op[node, :u_i] * dist_vec[1][2]])]
    elseif length(idx) == 3
        dist_vec_node = [sort([op[node, :u_r] * dist_vec[1][1], op[node, :u_r] * dist_vec[1][2]]), sort([op[node, :u_i] * dist_vec[1][1], op[node, :u_i] * dist_vec[1][2]]), dist_vec[2]]
    elseif length(idx) == 4
        dist_vec_node = [sort([op[node, :u_r] * dist_vec[1][1], op[node, :u_r] * dist_vec[1][2]]), sort([op[node, :u_i] * dist_vec[1][1], op[node, :u_i] * dist_vec[1][2]]), dist_vec[2], dist_vec[3]]
    else
        error("Wrong number of indices")
    end 
    return dist_vec_node
end

"""
    surv_vol(pg::PowerGrid, sol::ODESolution, u_r_idx::Array{Int,1}, u_i_idx::Array{Int,1}, dist_vec::Array{Float64,1})

Check it voltage is out of bounds for more than t_max seconds.
"""
function surv_vol(pg::PowerGrid, sol::ODESolution, op, u_r_idx::Array{Int,1}, u_i_idx::Array{Int,1})
    t = sol.t
    # The voltage of the Slack can not be perturbed_nodes
    # When the Slack has no imaginary part this could lead to problems
    # The voltage might shift from e-32 to e-21 and the state is wrongly counted
    for n in 1:length(pg.nodes)
        if typeof(pg.nodes[n]) == SlackAlgebraic
            continue
        end

        # difference between time series and operationpoint
        diff_v = abs.(1.0 .- (abs.(sol[u_r_idx[n], :] .+ 1im * sol[u_i_idx[n], :]) ./ op[n, :v]))
        a_v = findall(map(y -> y > 0.1, diff_v))  # Find all points in time series where they differ more than allowed
        if a_v != Int64[]
            if length(a_v) > 1
                cons_v = group_consecutive(a_v) # group consecutive points in time together
                deleteat!(cons_v, findall(length.(cons_v) .== 1))
                if length(cons_v) > 0
                    Δt_r = sum(map(k -> t[cons_v[k][end]] - t[cons_v[k][1]], 1:length(cons_v)))
                    if Δt_r > 5.0 # EN 50160 requirements
                        return 0.0
                    end
                end
            end
        end
    end    
    return 1.0
end

"""
    snbs_surv(pg::PowerGrid, rpg::ODEFunction, op, ω_idx, sim_time::Float64, numtries, dist_vec; tol_bs::Float64, tol_surv::Float64, tau_max::Float64)

A single trial for the Basin Stability and the Survivability. Not intended for direct usage. Use parallel_snbs_surv() instead.
"""
function snbs_surv(pg::PowerGrid, rpg::ODEFunction, g, op::State, ω_idx, ur_idx, ui_idx, indices, dist_vec_nodes, sim_time::Float64, numtries, dist_vec, tol_bs::Float64, tau_max::Float64; tol = false, write_frand = false)
    if write_frand
        file_name = "Frands_" * "Nodes_" * string(length(pg.nodes)) * ".csv"
        output = joinpath(@__DIR__, "../data/", file_name)
    end
    num_nodes = length(pg.nodes)

    # allocate return variables
    snbs = zeros(num_nodes)
    snbs_ω = zeros(num_nodes)
    snbs_u = zeros(num_nodes)

    surv = zeros(num_nodes)
    surv_ω = zeros(num_nodes)
    surv_u = zeros(num_nodes)

    infeasible = zeros(num_nodes)

    # run over all nodes
    for node in 1:num_nodes
        idx = indices[node] # index of the variables we want to perturb
        dist_vec_node = dist_vec_nodes[node]

        # drawing a random perturbation and simulate it
        rand_perturbation = random_perturbation(idx, rpg, op, dist_vec_node, tau_max)
        sol = simulate_perturbation(rpg, rand_perturbation, (0.0, sim_time), tol = tol)

        # only successful simulations are evaluated
        if sol.retcode == :Success
            final_state = sol.u[end]

            max_final_frequency = maximum(abs.(final_state[ω_idx])) # accessing the final state and get the maximal frequency
            final_diff_v = abs.(op[:, :v] .- (abs.(final_state[ur_idx] .+ 1im * final_state[ui_idx]) ./ op[:, :v]))

            max_frequency = maximum(vcat(sol[ω_idx,:]...)) # access the whole time series and get the maximal frequency
            min_frequency = minimum(vcat(sol[ω_idx,:]...))

            if max_final_frequency < tol_bs
                snbs_ω[node] = 1.0 
            end
            
            if maximum(final_diff_v) < 0.1
                snbs_u[node] = 1.0
            end

            if snbs_ω[node] == 1.0 && snbs_u[node] == 1.0
                snbs[node] = 1.0
            end 

            if max_frequency < dist_vec[2][2] && min_frequency > dist_vec[2][1]
                surv_ω[node] = 1.0
            end

            surv_u[node] = surv_vol(pg, sol, op, ur_idx, ui_idx)

            if surv_ω[node] == 1.0 && surv_u[node] == 1.0
                surv[node] = 1.0
            end

            # save frand when only the omega condition is violated
            if write_frand
                if snbs_ω[node] == 0.0 && snbs_u[node] == 1.0
                    df = DataFrame(F_RAND = [rand_perturbation], NODE = node)
                    if isfile(output)
                        CSV.write(output, df, append = true)
                    else
                        CSV.write(output, df, append = false)
                    end
                end 
            end
            
        elseif sol.retcode == :Unstable
            g_end = sum(abs.(g(sol.u[end])))
            if isnan(g_end)
                if length(sol.t) > 1
                    if sum(abs.(g(sol.u[end - 1]))) > 0.0001 
                        infeasible[node] = 1.0
                    end
                end
            else
                if g_end > 0.0001 
                    infeasible[node] = 1.0
                end
            end
        end
    end
    return snbs, snbs_u, snbs_ω, surv, surv_u, surv_ω, infeasible
end

"""
    parallel_snbs_surv(pg::PowerGrid, sim_time::Float64, sample_size::Int, dist_vec, tol_bs::Float64, tol_surv::Float64, tau_max::Float64, method::String)

Calculates the Basin Stability and the Survivability. Uses pmap and snbs_surv() for the parallelization.
"""
function parallel_snbs_surv(pg::PowerGrid, op::State, sim_time::Float64, sample_size::Int, 
                            dist_vec, tol_bs::Float64,
                            tau_max::Float64, method::String, save::Bool;
                            tol = false, write_frand = false)
    file_name =  method * "_Nodes_" * string(length(pg.nodes)) * "_SampleSize_" * string(sample_size) * "_TolBS_" * string(tol_bs) * "_SimTime_" * string(sim_time) * "_TauMax_" * string(tau_max) * "_Dist_" * string(dist_vec) * ".txt"
    rpg, g, ω_idx, ur_idx, ui_idx, indices, dist_vec_nodes = preliminary_snbs_surv(pg, method, op, dist_vec) # preliminary calculations we only need to do once
    
    result = pmap(sample -> snbs_surv(pg, rpg, g, op, ω_idx, ur_idx, ui_idx, indices, dist_vec_nodes, sim_time, sample, dist_vec, tol_bs, tau_max; tol = tol, write_frand = write_frand), 1:sample_size)

    snbs = mean(map(x -> x[1], result))
    snbs_u = mean(map(x -> x[2], result))
    snbs_ω = mean(map(x -> x[3], result))
    surv = mean(map(x -> x[4], result))
    surv_u = mean(map(x -> x[5], result))
    surv_ω = mean(map(x -> x[6], result))
    infeasible = mean(map(x -> x[7], result))

    if save
        path = joinpath(@__DIR__, "../data/")

        open(path * "SNBS_SURV_" * file_name, "w") do io
            writedlm(io, [snbs snbs_u snbs_ω surv surv_u surv_ω infeasible])
        end
        return nothing
    else
        return snbs, snbs_u, snbs_ω, surv, surv_u, surv_ω, infeasible
    end
end

"""
    preliminary_snbs_surv(pg::PowerGrid)

Calculates some preliminary thing which stay the same during all runs.
"""
function preliminary_snbs_surv(pg::PowerGrid, method::String, op, dist_vec)
    rpg = rhs(pg) # right hand side function
    g = constraint_equations(rpg) # constraint equations

    ω_idx = findall(map(variable -> occursin("ω", variable), string.(rpg.syms))) # indexes of ω variables, needs to work for Arrays and Dicts
    ur_idx = findall(map(variable -> occursin("u_r", variable), string.(rpg.syms)))
    ui_idx = findall(map(variable -> occursin("u_i", variable), string.(rpg.syms)))
    indices = map(n -> pd_node_idx(pg, rpg, n, method), 1:length(pg.nodes))
    dist_vec_nodes = map(n -> adapted_dist_vec(indices[n], n, op, dist_vec), 1:length(pg.nodes))

    return rpg, g, ω_idx, ur_idx, ui_idx, indices, dist_vec_nodes
end