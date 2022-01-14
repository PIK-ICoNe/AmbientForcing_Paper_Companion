
gr()
"""
    plot_res(result::PowerGridSolution)

Plots the result of a PowerGridSolution.
"""
function plot_res(result::PowerGridSolution)
    powergrid = result.powergrid
    ω_indices = findall(n -> :ω ∈ symbolsof(n), powergrid.nodes);

    n = length(powergrid.nodes)

    CList = reshape( range(colorant"coral1", stop=colorant"steelblue",length=n), 1, n );

    ω_labels = reshape([latexstring(string(raw"\omega", "_{$i}")) for i=ω_indices], (1, length(ω_indices)))
    #p_labels = reshape([latexstring(string(raw"p", "_{$i}")) for i=1:length(powergrid.nodes)], (1, length(powergrid.nodes)))
    v_labels = reshape([latexstring(string(raw"v", "_{$i}")) for i=1:length(powergrid.nodes)], (1, length(powergrid.nodes)))

    pl_ω = Plots.plot(result, ω_indices, :ω, linecolor = CList, legend = false, ylabel=L"\omega \left[rad/s\right]", label=ω_labels)
    #pl_p = Plots.plot(result, :, :p, linecolor = CList, legend = false, ylabel=L"p [p.u.]", label=p_labels)
    pl_v = Plots.plot(result, :, :v, linecolor = CList, legend = false, ylabel=L"v [p.u.]", label=v_labels)

    plt = Plots.plot(
        pl_ω,pl_v;
        layout=(2,1),
        size = (500, 500),
        lw=3,
        xaxis = (L"t[s]")
    )
    display(plt)
    return plt
end

"""
    plot_vs(pg::PowerGrid, data1, data2, err1, err2, xlab, ylab, my_color)
"""
function plot_vs(pg::PowerGrid, data1, data2, err1, err2, xlab, ylab, color1, color2, legend, lims)
    # Differentiate between dynamical and Algebraic
    cons_idx, dyn_idx, constraint_vec = constraint_vector(pg)
    Plots.plot(data1[cons_idx],
            data2[cons_idx],
            xerror = err1[cons_idx],
            yerror = err2[cons_idx],
            xlims = lims,
            ylims = lims,
            seriestype = :scatter,
            xaxis = xlab,
            yaxis = ylab,
            linecolor = color1,
            opacity = 0.8,
            marker = (:hexagon, 10, color1),
            label = "Algebraic",
            show = true
    )

    Plots.plot!(data1[dyn_idx],
            data2[dyn_idx],
            xerror = err1[dyn_idx],
            yerror = err2[dyn_idx],
            seriestype = :scatter,
            opacity = 0.8,
            linecolor = color2,
            marker = (:circle, 10, color2),
            foreground_color_legend = nothing,
            label = "Differential",
            legend = legend,
            show = true
    )

    plt = Plots.plot!([0,1],
            [0,1],
            line = :dash,
            label = "",
            color = colorant"gray0"
    )
    return plt
end

function plot_vs_only_dyn(pg::PowerGrid, data1, data2, err1, err2, xlab, ylab, my_color, lims)
    cons_idx, dyn_idx, constraint_vec = constraint_vector(pg)
    Plots.plot(data1[dyn_idx],
            data2[dyn_idx],
            xerror = err1[dyn_idx],
            yerror = err2[dyn_idx],
            seriestype = :scatter,
            xlims = lims,
            ylims = lims,
            xaxis = xlab,
            yaxis = ylab,
            opacity = 0.8,
            linecolor = my_color,
            #series_annotations = ([string(x) for x in dyn_idx]),
            marker = (:circle, 10, my_color),
            foreground_color_legend = nothing,
            legend = false,
            show = true
    )

    plt = Plots.plot!([0,1],
            [0,1],
            line = :dash,
            label = "",
            color = colorant"gray0"
    )
    return plt
end

"""
    spread_violin_plot(distances_sorted, spread_sorted, color1, color2, dist_measure)
"""
function spread_violin_plot(distances_sorted, spread_sorted, color1, color2, dist_measure)
    StatsPlots.violin(string.(Int.(distances_sorted)), spread_sorted, labelfontsize = 15, 
                      color = color1, xlabel = dist_measure,
                      ylabel = L"s(k)", linecolor = color1 , outliers = true,
                      label = "Violin Plot")

    plt = StatsPlots.boxplot!(string.(Int.(distances_sorted)), spread_sorted, labelfontsize = 15, 
                        color = color2, xlabel = dist_measure,
                        ylabel = L"s(k)", linewidth = 2, fillalpha = 0.5, linecolor = :black, 
                        outliers = false, label = "Box Plot")
    display(plt)
    return plt
end

"""
    my_corrplot(df::DataFrame)
"""
function my_corrplot(df::DataFrame, data1::String, data2::String, legendpos)
    Classes = ["Bulk", "Root", "Inner Tree Node", "Proper Leave", "Sparse Sprout", "Dense Sprout"]
    df_classes = []

    # Make a filtered DataFrame for each Node Class
    for i in Classes
        append!(df_classes, [filter("NodeClass" => x -> x == i, df)])
    end

    # setting the layout for the 3 different plots
    layout = @layout [a            _
                    b{0.8w,0.8h} c]

    plt = Plots.plot(layout = layout, link = :both, size = (500, 500), margin = -10Plots.px)
    Plots.plot!([0,1], [0,1], line = :dash, subplot = 2,label = "", foreground_color_legend = nothing,color = colorant"gray42", legend = legendpos)
    color_vec = distinguishable_colors(6, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)

    # Add a separate histogram and scatter plot for each class in a different color
    for j in 1:length(df_classes)
        Plots.scatter!(df_classes[j][!,data1], legendfontsize = 8, df_classes[j][!,data2], label = Classes[j], left_margin = 1mm, xlabel = L"\beta", ylabel = L"\sigma", subplot = 2, xlims = [0,1], ylims = [0,1])
        Plots.histogram!([df_classes[j][!,data1] df_classes[j][!,data2]], bins = 20, subplot = [1 3], orientation = [:v :h], framestyle = :none, linewidth = 0, legend = false)
    end
    display(plt)
end

"""
    function my_graph_plot(pg::PowerGrid, df::DataFrame, pg_idx::Int, label_nodes = [])

Using GraphMakie to plot some nice powergrids. Markersize is with respect to the SNBS. 
"""
function my_graph_plot(pg::PowerGrid, df::DataFrame, pg_idx::Int, lable_nodes = [])
    df_pg = df[df.PG_IDX .== pg_idx, :] # use only the subset of the data frame

    num_nodes = length(pg.nodes)
    pa = palette(:tab10)
    node_color = fill(colorant"grey", num_nodes)
    node_marker = fill(:utriangle, num_nodes)

    # Give each Class a specific Color
    for i in 1:num_nodes
        class = df_pg[i, :NodeClass]
        if class == "Bulk"
            node_color[i] = pa[1]
        elseif class == "Root"
            node_color[i] = pa[2]
        elseif class == "Inner Tree Node"
            node_color[i] = pa[3]
        elseif class == "Proper Leave"
            node_color[i] = pa[4]
        elseif class == "Sparse Sprout"
            node_color[i] = pa[5]
        elseif class == "Dense Sprout"
            node_color[i] = pa[6]
        end
        if df_pg[i, :Constraint] == 0.0
            node_marker[i] = :circle
        end
    end

    if lable_nodes != []
        node_lable = fill("", num_nodes)
        if length(lable_nodes) == 1
            node_lable[lable_nodes] = string(lable_nodes)
        else    
            map(n -> node_lable[n] = string(n), lable_nodes)
        end
        f, ax, p = graphplot(pg.graph, node_marker = node_marker, nlabels = node_lable, node_color = node_color, node_size = df_pg[!, :SNBS] .* 20)
    else
        f, ax, p = graphplot(pg.graph, node_marker = node_marker, node_color = node_color, node_size = df_pg[!, :SNBS] .* 20)
    end
    #hidedecorations!(ax); hidespines!(ax)
    #ax.aspect = DataAspect()
    return f
end 

"""
    density_class_plot(df::DataFrame, data::String; legend_pos = :topleft)
"""
function density_class_plot(df::DataFrame, data::String, xaxis; legend_pos = :topleft)
    pa = palette(:tab10)
    StatsPlots.density(df[df.NodeClass .== "Bulk", data], color = pa[1], lw = 3, xaxis = xaxis, xlims = [0, 1], label = "Bulk", legendfontsize = 11, legend = legend_pos)
    StatsPlots.density!(df[df.NodeClass .== "Root", data], color = pa[2], lw = 3, label = "Root")
    StatsPlots.density!(df[df.NodeClass .== "Inner Tree Node", data], color = pa[3], lw = 3, label = "Inner Tree Node")
    StatsPlots.density!(df[df.NodeClass .== "Proper Leave", data], color = pa[4], lw = 3, label = "Proper Leave")
    StatsPlots.density!(df[df.NodeClass .== "Sparse Sprout", data], color = pa[5], lw = 3, label = "Sparse Sprout")
    StatsPlots.density!(df[df.NodeClass .== "Dense Sprout", data], color = pa[6], lw = 3, label = "Dense Sprout")
end

"""
    density_degree_plot(df::DataFrame, data::String, xaxis; xlims = [0,1], color = ColorSchemes.berlin, legend_pos = :topleft)
"""
function density_degree_plot(df::DataFrame, data::String, xaxis, yaxis; xlims = [0,1], color = ColorSchemes.berlin, legend_pos = :topleft)
    max_degree = maximum(df[!, "degree"])
    my_range = 0:0.01:1

    nodes_per_degree = map(i -> nrow(filter("degree" => x -> x == i, df)), 1:max_degree)
    last_degree = findlast(nodes_per_degree .> 10)

    mi = range(0,1, length = last_degree)
    pa = color[mi]

    dens = kde(df[df.degree .== 1, data])
    dens = Spline1D(dens.x, dens.density)

    plt = plot(collect(my_range), dens(my_range) ./ sum(dens(my_range)), label = "1", color = pa[1], lw = 3, xaxis = xaxis, yaxis = yaxis, xlims = xlims, legendfontsize = 11, legend = legend_pos)

    for i in 2:last_degree
        dens = kde(df[df.degree .== i, data])
        dens = Spline1D(dens.x, dens.density)

        plt = plot!(collect(my_range), dens(my_range) ./ sum(dens(my_range)), label = string(i), color = pa[i], lw = 3)
    end
    display(plt)
end


"""
    density_class_plot!(df::DataFrame, data::String; legend_pos = :topleft)
"""
function density_degree_plot!(df::DataFrame, data::String, xaxis; xlims = [0,1], color = ColorSchemes.berlin, legend_pos = :topleft)
    max_degree = maximum(df[!, "degree"])

    nodes_per_degree = map(i -> nrow(filter("degree" => x -> x == i, df)), 1:max_degree)
    last_degree = findlast(nodes_per_degree .> 10)

    mi = range(0,1, length = last_degree)
    pa = color[mi]


    plt = StatsPlots.density!(df[df.degree .== 1, data], label = "", color = pa[1], lw = 3, xaxis = xaxis, xlims = xlims, legendfontsize = 11, legend = legend_pos)

    for i in 2:last_degree
        StatsPlots.density!(df[df.degree .== i, data], label = "", color = pa[i], lw = 3)
    end

    display(plt)
end

"""
    plot_fault_simulation(pg::PowerGrid, rpg::ODEFunction, op::State, sim_time::Float64, pert_vec, node::Int)
"""
function plot_fault_simulation(pg::PowerGrid, rpg::ODEFunction, op::State, sim_time::Float64, pert_vec, node::Int)
    pert_vec = [pert_vec[1] * op[node, :u_r], pert_vec[2] * op[node, :u_i], pert_vec[3]]
    z = set_perturbation(pg, rpg, node, op, pert_vec)
    sol = simulate_perturbation(rpg, z, (0.0, sim_time))
    plot_res(PowerGridSolution(sol, pg))
end