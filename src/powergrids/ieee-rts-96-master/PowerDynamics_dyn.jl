function get_ieee_96()
    ### Load Network Structure and Parameters from Data
    path = string(@__DIR__) * "/Data/"

    BusData = CSV.read(joinpath(path,"Bus.csv"), DataFrame)
    LineData = CSV.read(joinpath(path,"Line.csv"), DataFrame)
    LineData[!, :source] = LineData.source .|> Int
    LineData[!, :dest] = LineData.dest .|> Int
    GeneratorData = CSV.read(joinpath(path,"Generator.csv"), DataFrame)
    LoadData = CSV.read(joinpath(path,"Load.csv"), DataFrame)
    FlowData = CSV.read(joinpath(path,"Flow.csv"), DataFrame)

    N = nrow(BusData)
    L = nrow(LineData)
    G = nrow(GeneratorData)

    g = SimpleGraph(N)
    for i = 1:L
        add_edge!(g, Int64(LineData[i, :source]), Int64(LineData[i, :dest]))
    end

    node_df = outerjoin(BusData, GeneratorData, LoadData, FlowData, on=:ID, makeunique = true)

    slack_idx = argmax(skipmissing(node_df.P_Gen))

    nodes = []
    for n in eachrow(node_df)
        if n.Number == slack_idx
            push!(nodes, SlackAlgebraic(U = 1.0))
        else
            if n.P_Gen |> ismissing
                push!(nodes, PVAlgebraic(P = -n.P_Load, V = 1.0))
            else
                push!(nodes, PVAlgebraic(P = n.P_Gen - n.P_Load, V = 1.0))
            end
        end
    end

    lines = []
    for l in eachrow(LineData)
        if node_df[l.source, :Base_V] == node_df[l.dest, :Base_V] # Normal line
            push!(lines, StaticLine(; from = l.source, to = l.dest, Y = inv(complex(l.r, l.x))))
        else
            t_ratio = node_df[l.dest, :Base_V] / node_df[l.source, :Base_V]
            push!(lines, Transformer(; from = l.source, to = l.dest, y = inv(complex(l.r, l.x)), t_ratio = t_ratio) )
        end
    end
    pg_cons = PowerGrid(nodes, lines)
    op_cons = find_operationpoint(pg_cons)

    nodes = []
    for i in 1:length(pg_cons.nodes)
        n = node_df[i, :]
        if n.Number == slack_idx
            push!(nodes, SlackAlgebraic(U=1.0)) 
        else
            if n.P_Gen |> ismissing
                push!(nodes, PQAlgebraic(P = op_cons[i, :p], Q = op_cons[i, :q]))
            else
                push!(nodes, SchifferApprox(τ_P = n.Inertia, τ_Q = 8.0, K_P = 5, K_Q = 0.1, V_r = 1.0, P = op_cons[i, :p], Q = op_cons[i,:q], Y_n = 0))
            end
        end
    end
    PowerGrid(nodes, lines)
end