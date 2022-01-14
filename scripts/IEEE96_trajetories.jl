using Pkg
Pkg.activate(joinpath(@__DIR__, "../"))
Pkg.instantiate()
using BS_DAE

pg_ieee96 = get_ieee_96()

op = find_operationpoint(pg_ieee96)
rpg_96 = rhs(pg_ieee96)

sim_time = 200.0
sol_op = simulate_perturbation(rpg_96, op.vec, (0.0, sim_time))
g = constraint_equations(rhs(pg_ieee96))
λ, stable = check_eigenvalues(pg_ieee96, op)

node = 2
pert_vec_up = [abs(op[node, :u_r]) * -0, abs(op[node, :u_i]) * -0, 3.0]
z_up = set_perturbation(pg_ieee96, rpg_96, node, op, pert_vec_up)
sol_up = simulate_perturbation(rpg_96, z_up, (0.0, sim_time))

vect = map(x -> log(sum(abs.(g(x)))), sol_up.u)
vec2 = map(x -> log(norm(x)), sol_up.u)

g(sol_up.u[end - 1])
plot(sol_up.t, vect, legend = false)



#nodes_up = perturbed_nodes(rpg_96, op, sol_up)
plot_res(PowerGridSolution(sol_up, pg_ieee96))
png(joinpath(@__DIR__,"../plots/ieee96/example_trajectory_up_simtime" * string(sim_time) * "_node_" * string(node)))

plot_res_v(PowerGridSolution(sol_up, pg_ieee96), 1:73 ,colorant"coral1", colorant"steelblue")
plot_res_v!(PowerGridSolution(sol_op, pg_ieee96), nodes_up ,colorant"coral1", colorant"steelblue")
png(joinpath(@__DIR__,"../plots/ieee96/trajectory/example_trajectory_up_vol_simtime" * string(sim_time) * "_node_" * string(node)))


#=
pert_vec_down = [abs(op[node, :u_r]) * -0.1, abs(op[node, :u_i]) * -0.1, -2.0]
z_down = set_perturbation(pg_ieee96, rpg_96, node, op, pert_vec_down)
sol_down = simulate_perturbation(rpg_96, z_down, (0.0, sim_time))
nodes_down = perturbed_nodes(rpg_96, op, sol_down)

plot_res(PowerGridSolution(sol_down, pg_ieee96), nodes_down, colorant"coral1", colorant"steelblue")
png(joinpath(@__DIR__,"../plots/ieee96/example_trajectory_down_simtime" * string(sim_time) * "_node_" * string(node)))

plot_res_v(PowerGridSolution(sol_down, pg_ieee96), nodes_down, colorant"coral1", colorant"steelblue")
plot_res_v!(PowerGridSolution(sol_op, pg_ieee96), nodes_down, colorant"coral1", colorant"steelblue")
png(joinpath(@__DIR__,"../plots/ieee96/trajectory/example_trajectory_down_vol_simtime" * string(sim_time) * "_node_" * string(node)))

=#