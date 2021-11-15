using Test, FinalProjectPedroFerraz
using Graphs, SimpleWeightedGraphs, StatsBase, Random, DataFrames, XLSX, BenchmarkTools, Plots, Query

dataset_path = (@__DIR__) * "/data/goeke-2018/"
file_name = "rc204C6" * ".txt"
function read_file(file_name, distance_scaling = 1.0)
    file = readchomp(dataset_path * file_name)
    lines = split(file, r"\n")
    rows = split.(lines[1:end-6], r" +")
    cols = [[rows[i][j] for i in 2:length(rows)] for j in 1:length(rows[1])]
    df = DataFrame(cols[1:end], rows[1])
    # df = DataFrame(x = parse.(Float64, df.x), y = parse.(Float64, df.y), 
    #                 readytime = parse.(Float64, df.ReadyTime))

    origin_df = @from i in df begin
        @where i.Type == "d"
        @select {i.x, i.y}
        @collect DataFrame
    end
    x_coord = parse.(Float64, origin_df.x)[1]
    y_coord = parse.(Float64, origin_df.y)[1]

    requests_df = @from i in df begin
        @where i.Type == "cd" || i.Type == "cp"
        @select {i.x, i.y, i.ReadyTime}
        @collect DataFrame
    end
    x_coords = parse.(Float64, requests_df.x)
    y_coords = parse.(Float64, requests_df.y)
    release_times = parse.(Float64, requests_df.ReadyTime)
    
    points = hcat([x_coord, y_coord], vcat([x_coords', y_coords']...))
    points ./= distance_scaling
    return points, release_times
end

function build_instance(file_name, distance_scaling = 1.0)
    points, release_times = read_file(file_name, distance_scaling)
    g, dists = euclidean_graph(points, p=2)
    srcs = Vector{Int}()
    dsts = Vector{Int}()
    wgts = Vector{Float64}()
    for e in dists
        push!(srcs, e.first.src)
        push!(dsts, e.first.dst)
        push!(wgts, e.second)
    end
    g = SimpleWeightedGraph(srcs, dsts, wgts);
    
    N_req = length(release_times)
    request_vertices = 2:nv(g)
    requests = map(i -> DataModel.Request(release_times[i], request_vertices[i]), 1:N_req)

    return g, requests
end

function run_online_tests(instance_names)
    main_df = DataFrame()
    for (index, instance_name) in enumerate(instance_names)
        # TODO: colocar uma "barrinha de progresso"
        dfs = []
        g, requests = build_instance(instance_name * ".txt")
        should_set_header = true
        capacities = filter(c -> c < nv(g), [1, 2, 3, 5])
        for capacity in capacities
            local cost_ignore, end_t_ignore, route_ignore
            local cost_return, end_t_return, route_return
            local cost_naive, end_t_naive, route_naive
            local cost_naive_return, end_t_naive_return, route_naive_return
            local cost_compute_return, end_t_compute_return, route_compute_return
            local exec_time_ignore, exec_time_return, exec_time_naive, exec_time_naive_return, exec_time_compute_return
            skip = false

            try
                # TODO: permitir que apenas o Wait and Return tenha resultados, por exemplo
                ((cost_naive, end_t_naive, route_naive), exec_time_naive) = @timed NaiveIgnore.run(g, requests, capacity, 8, 25)
                ((cost_ignore, end_t_ignore, route_ignore), exec_time_ignore) = @timed WaitAndIgnore.run(g, requests, capacity, 8, 25)
                ((cost_return, end_t_return, route_return), exec_time_return) = @timed WaitAndReturn.run(g, requests, capacity, 8, 25)
                ((cost_naive_return, end_t_naive_return, route_naive_return), exec_time_naive_return) = @timed NaiveReturn.run(g, requests, capacity, 8, 25)
                ((cost_compute_return, end_t_compute_return, route_compute_return), exec_time_compute_return) = @timed ComputeReturn.run(g, requests, capacity, 8, 25)
            catch e
                println(e)
                skip = true
            end
            if skip
                continue
            end
            costs = [cost_ignore, cost_return, cost_naive, cost_naive_return, cost_compute_return]
            end_times = [end_t_ignore, end_t_return, end_t_naive, end_t_naive_return, end_t_compute_return]
            exec_times = [exec_time_ignore, exec_time_return, exec_time_naive, exec_time_naive_return, exec_time_compute_return]
            routes = repr.([route_ignore, route_return, route_naive, route_naive_return, route_compute_return])

            instance_element = missing
            if should_set_header
                instance_element = "$instance_name: $(nv(g)) nós, $(length(requests)) pedidos"
                should_set_header = false
            end
            instance_array = [instance_element, repeat([missing], length(costs)-1)...]
            capacity_array = [capacity, repeat([missing], length(costs)-1)...]
            algorithms = ["Wait and Ignore", "Wait and Return", "Naive Ignore", "Naive Return", "Compute Return"]

            df = DataFrame("Instância" => instance_array,
                            "Capacidade" => capacity_array,
                            "Algoritmo" => algorithms,
                            "Custo" => costs,
                            "Tempo de término" => end_times,
                            "Rota encontrada" => routes,
                            "Tempo de execução" => exec_times
                        )
            push!(dfs, df)
        end

        if index == 1
            main_df = dfs[1]
            for df in dfs[2:end]
                append!(main_df, df)
            end
        else
            for df in dfs
                append!(main_df, df)
            end
        end
    end
    path = "test/results/"
    mkpath(path)
    XLSX.writetable(path * "results_large.xlsx", collect(DataFrames.eachcol(main_df)), DataFrames.names(main_df), overwrite=true)
    main_df
end

function plot_instance(instance_name, distance_scaling = 1.0)
    points, release_times = read_file(instance_name * ".txt", distance_scaling)
    plot = Plots.scatter(points[1, :], points[2, :], markersize=6, legend=false)
    ydiff = maximum(points[2, :]) - minimum(points[2, :])
    annotate!(points[1, 1], points[2, 1], text(1, 4))
    for i in 2:size(points, 2)
        annotate!(points[1, i], points[2, i], text(i, 4))
        annotate!(points[1, i], points[2, i] + 0.03 * ydiff, text(Int64(release_times[i-1]), 3))
    end
    path = "test/images/"
    mkpath(path)
    Plots.savefig(plot, path * "$instance_name.pdf")
end

instances = readdir(dataset_path)[1:end-1]
instance_names = first.(split.(instances, "."))
@time df = run_online_tests(instance_names)

plot_instance("rc204C6")
points, release_times = read_file("rc204C6" * ".txt")
size(points, 2)

plot_instance("c101C6")

g, requests = build_instance("lrc101" * ".txt")
MilpOfflineAlgorithm.run(g, requests, 5, full_search=false, time_limit=10.0)
ComputeReturn.run(g, requests, 5, 8, 25)
NaiveIgnore.run(g, requests, 5, 8)
