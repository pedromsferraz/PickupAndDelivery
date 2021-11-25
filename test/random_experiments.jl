using Test, FinalProjectPedroFerraz
using Graphs, SimpleWeightedGraphs, StatsBase, Distributions, DataFrames, XLSX, BenchmarkTools, Plots, Plots.PlotMeasures

# Generate an euclidean graph with N request uniformly distributed in the box [0, L]^2
# with requests distributed ~ Exp(β)
function build_random_instance(N, L, β)
    # Create an euclidean graph 
    g, dists, points = euclidean_graph(N+1, 2, p=2, L=L)
    srcs = Vector{Int}()
    dsts = Vector{Int}()
    weights = Vector{Float64}()
    for e in dists
        push!(srcs, e.first.src)
        push!(dsts, e.first.dst)
        push!(weights, e.second)
    end
    g = SimpleWeightedGraph(srcs, dsts, weights);

    # Create requests
    t = Exponential(β)
    wait_times = rand(t, N)
    release_times = [sum(wait_times[1:n]) for n in 1:N]
    request_vertices = sample(2:nv(g), N, replace=false)
    requests = map(i -> DataModel.Request(release_times[i], request_vertices[i]), 1:N)
    
    return g, requests, points, release_times
end

function plot_instance(points, requests)
    plot = Plots.scatter(points[1, :], points[2, :], markersize=6, legend=false)
    ydiff = maximum(points[2, :]) - minimum(points[2, :])
    annotate!(points[1, 1], points[2, 1], text(1, 4))
    for i in 2:size(points, 2)
        annotate!(points[1, i], points[2, i], text(requests[i-1].destination, 4))
        annotate!(points[1, i], points[2, i] + 0.03 * ydiff, text(floor(requests[i-1].release_time), 5))
    end
    path = "test/images/"
    mkpath(path)
    Plots.savefig(plot, path * "Random.pdf")
end

function run_competitive_ratio_random_tests(k, N, L, Λ)
    main_df = DataFrame()
    for (index, β) in enumerate(Λ)
        println("β = $β")
        # TODO: colocar uma "barrinha de progresso"
        dfs = []
        should_set_header = true
        instances = [build_random_instance(N, L, β) for i in 1:k]
        g, _ = instances[1]
        capacities = filter(c -> c < nv(g), [1, 2, 3, 5])
        for capacity in capacities
            println("c = $capacity")
            costs = [0., 0., 0., 0., 0., 0.]
            end_times = [0., 0., 0., 0., 0., 0.]
            exec_times = [0., 0., 0., 0., 0., 0.]
            competitive_ratios = [0., 0., 0., 0., 0., 0.]
            number_of_runs = 0
            
            for i in 1:k
                g, requests = instances[i]
                local cost_ignore, end_t_ignore, route_ignore
                local cost_return, end_t_return, route_return
                local cost_naive, end_t_naive, route_naive
                local cost_naive_return, end_t_naive_return, route_naive_return
                local cost_compute_return, end_t_compute_return, route_compute_return
                local cost_offline, end_t_offline, route_offline
                local exec_time_ignore, exec_time_return, exec_time_naive
                local exec_time_naive_return, exec_time_compute_return, exec_time_offline
                skip = false
                
                while !skip
                    try
                        ((cost_naive, end_t_naive, route_naive), exec_time_naive) = @timed NaiveIgnore.run(g, requests, capacity, 8)
                        ((cost_ignore, end_t_ignore, route_ignore), exec_time_ignore) = @timed WaitAndIgnore.run(g, requests, capacity, 8)
                        ((cost_return, end_t_return, route_return), exec_time_return) = @timed WaitAndReturn.run(g, requests, capacity, 8)
                        ((cost_naive_return, end_t_naive_return, route_naive_return), exec_time_naive_return) = @timed NaiveReturn.run(g, requests, capacity, 8)
                        ((cost_compute_return, end_t_compute_return, route_compute_return), exec_time_compute_return) = @timed ComputeReturn.run(g, requests, capacity, 8)
                        ((cost_offline, end_t_offline, route_offline), exec_time_offline) = @timed OfflineAlgorithm.run_multithreaded(g, requests, capacity)
                        break
                    catch e
                        println(e.msg)
                        skip = true
                    end
                end
                if skip
                    continue
                end

                number_of_runs += 1
                current_costs = [cost_ignore, cost_return, cost_naive, cost_naive_return, cost_compute_return, cost_offline]
                costs += current_costs
                end_times += [end_t_ignore, end_t_return, end_t_naive, end_t_naive_return, end_t_compute_return, end_t_offline]
                exec_times += [exec_time_ignore, exec_time_return, exec_time_naive, exec_time_naive_return, exec_time_compute_return, exec_time_offline]
                competitive_ratios += current_costs ./= cost_offline
            end

            if number_of_runs == 0
                continue
            end
            costs /= number_of_runs
            end_times /= number_of_runs
            exec_times /= number_of_runs
            competitive_ratios /= number_of_runs
    
            instance_element = missing
            if should_set_header
                instance_element = "N=$N, L=$L, β=$β"
                should_set_header = false
            end
            instance_array = [instance_element, repeat([missing], length(costs)-1)...]
            capacity_array = [capacity, repeat([missing], length(costs)-1)...]
            number_of_runs_array = [number_of_runs, repeat([missing], length(costs)-1)...]
            algorithms = ["Wait and Ignore", "Wait and Return", "Naive Ignore", "Naive Return", "Compute Return", "Offline"]
    
            df = DataFrame("Instância" => instance_array,
                            "Número de execuções (k)" => number_of_runs_array,                
                            "Capacidade" => capacity_array,
                            "Algoritmo" => algorithms,
                            "Competitive Ratio dentre k execuções" => competitive_ratios,
                            "Custo médio dentre k execuções" => costs,
                            "Tempo de término médio dentre k execuções" => end_times,
                            "Tempo de execução dentre k execuções" => exec_times,
                        )
            push!(dfs, df)
        end
        if length(dfs) == 0
            continue
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
    XLSX.writetable(path * "competitive_ratio $N, $L, $Λ.xlsx", collect(DataFrames.eachcol(main_df)), DataFrames.names(main_df), overwrite=true)
    main_df
end

function run_online_random_tests(k, N, L, Λ)
    main_df = DataFrame()
    for (index, β) in enumerate(Λ)
        println("β = $β")
        # TODO: colocar uma "barrinha de progresso"
        dfs = []
        should_set_header = true
        instances = [build_random_instance(N, L, β) for i in 1:k]
        g, _ = instances[1]
        capacities = filter(c -> c < nv(g), [1, 2, 3, 5])
        for capacity in capacities
            println("c = $capacity")
            costs = [0., 0., 0., 0., 0.]
            end_times = [0., 0., 0., 0., 0.]
            exec_times = [0., 0., 0., 0., 0.]
            wins = [0, 0, 0, 0, 0]
            number_of_runs = 0
            
            for i in 1:k
                g, requests, points, release_times = instances[i]
                local cost_ignore, end_t_ignore, route_ignore
                local cost_return, end_t_return, route_return
                local cost_naive, end_t_naive, route_naive
                local cost_naive_return, end_t_naive_return, route_naive_return
                local cost_compute_return, end_t_compute_return, route_compute_return
                local exec_time_ignore, exec_time_return, exec_time_naive, exec_time_naive_return, exec_time_compute_return
                skip = false
                
                while !skip
                    try
                        # TODO: permitir que apenas o Wait and Return tenha resultados, por exemplo
                        ((cost_naive, end_t_naive, route_naive), exec_time_naive) = @timed NaiveIgnore.run(g, requests, capacity, 8, 25)
                        ((cost_ignore, end_t_ignore, route_ignore), exec_time_ignore) = @timed WaitAndIgnore.run(g, requests, capacity, 8, 25)
                        ((cost_return, end_t_return, route_return), exec_time_return) = @timed WaitAndReturn.run(g, requests, capacity, 8, 25)
                        ((cost_naive_return, end_t_naive_return, route_naive_return), exec_time_naive_return) = @timed NaiveReturn.run(g, requests, capacity, 8, 25)
                        ((cost_compute_return, end_t_compute_return, route_compute_return), exec_time_compute_return) = @timed ComputeReturn.run(g, requests, capacity, 8, 25)
                        break
                    catch e
                        println(e)
                        # skip = true
                    end
                end
                if skip
                    continue
                end
                
                # if !(cost_naive_return ≈ cost_compute_return) && cost_compute_return + 1.0 < cost_naive_return && !([0] in route_naive_return)
                #     println(cost_naive_return)
                #     println(cost_compute_return)
                #     return g, requests, points, release_times
                # end

                number_of_runs += 1
                current_costs = [cost_ignore, cost_return, cost_naive, cost_naive_return, cost_compute_return]
                costs += current_costs
                end_times += [end_t_ignore, end_t_return, end_t_naive, end_t_naive_return, end_t_compute_return]
                exec_times += [exec_time_ignore, exec_time_return, exec_time_naive, exec_time_naive_return, exec_time_compute_return]
                wins += [minimum(current_costs) ≈ current_costs[p] for p in 1:5]
            end

            if number_of_runs == 0
                continue
            end
            costs /= number_of_runs
            end_times /= number_of_runs
            exec_times /= number_of_runs
    
            instance_element = missing
            if should_set_header
                instance_element = "N=$N, L=$L, β=$β"
                should_set_header = false
            end
            instance_array = [instance_element, repeat([missing], length(costs)-1)...]
            capacity_array = [capacity, repeat([missing], length(costs)-1)...]
            number_of_runs_array = [number_of_runs, repeat([missing], length(costs)-1)...]
            algorithms = ["Wait and Ignore", "Wait and Return", "Naive Ignore", "Naive Return", "Compute Return"]
    
            df = DataFrame("Instância" => instance_array,
                            "Número de execuções (k)" => number_of_runs_array,                
                            "Capacidade" => capacity_array,
                            "Algoritmo" => algorithms,
                            "Melhor algoritmo dentre k execuções" => wins,
                            "Custo médio dentre k execuções" => costs,
                            "Tempo de término médio dentre k execuções" => end_times,
                            "Tempo de execução dentre k execuções" => exec_times,
                        )
            push!(dfs, df)
        end
        if length(dfs) == 0
            continue
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
    XLSX.writetable(path * "random_algorithms $N, $L, $Λ.xlsx", collect(DataFrames.eachcol(main_df)), DataFrames.names(main_df), overwrite=true)
    main_df
end

k = 100
N = 5
L = 500
Λ = [50, 100, 250, 500, 1000]
@time run_online_random_tests(k, N, L, Λ)

@time run_competitive_ratio_random_tests(k, N, L, Λ)




@time run_online_random_tests(10, 20, 500, [100, 1000])
plot_instance(points, requests)

g, requests, points, release_times = build_random_instance(20, 500, 100)
@time ComputeReturn.run(g, requests, 10, 8)

function benchmark_speedup(M)
    sum = 0
    for i in 1:M
        t1 = @elapsed OfflineAlgorithm.run(g, requests, 10)
        t2 = @elapsed OfflineAlgorithm.run_multithreaded(g, requests, 10)
        sum += t1 / t2
    end
    sum /= M
    sum
end
benchmark_speedup(100)

import Plots, XLSX
# RANDOM INSTANCE RESULTS
N = 5
data = XLSX.readdata("test/rand$N.xlsx", "Sheet1", "A1:H101")
capacity = data[2:5:end, 3]
data[2:25:end, 6]
wait_and_ignore = data[2:5:end, 6]
wait_and_return = data[3:5:end, 6]
naive_ignore = data[4:5:end, 6]
naive_return = data[5:5:end, 6]
compute_return = data[6:5:end, 6]

# Algorithms per β
# Avg Cost vs. Capacity for each β
function get_plots_for_index(index)
    β = [50, 100, 250, 500, 1000][index]
    c = capacity[1:4]
    plt = plot(c, [wait_and_ignore[1+4*k:4+4*k] for k in 0:4][index], label="Wait and Ignore", 
                title="β = $β", xlabel="capacidade", ylabel="custo médio", size = (900, 500), 
                xtickfontsize = 20, ytickfontsize = 20, titlefontsize = 20, labelfontsize = 20, 
                linewidth = 3, legendfontsize=20, left_margin = 5mm, bottom_margin=6mm)
    plot!(c, [wait_and_return[1+4*k:4+4*k] for k in 0:4][index], label="Wait and Return", linewidth = 3)
    plot!(c, [naive_ignore[1+4*k:4+4*k] for k in 0:4][index], label="Naive Ignore", linewidth = 3)
    plot!(c, [naive_return[1+4*k:4+4*k] for k in 0:4][index], label="Naive Return", linewidth = 3)
    plot!(c, [compute_return[1+4*k:4+4*k] for k in 0:4][index], label="Compute Return", linewidth = 3)
    return plt
end
for i in 1:5
    plt = get_plots_for_index(i)
    savefig(plt, "test/images/worst_caseN=$N,beta$([50, 100, 250, 500, 1000][i]).pdf")
end

# Algorithms per capacity
# Avg Cost vs. β for each capacity
function get_plots_for_index(index)
    c = capacity[1:4][index]
    β = [50, 100, 250, 500, 1000]
    plot(β, [wait_and_ignore[1+k:4:end] for k in 0:4][index], label="Wait and Ignore", 
            title="c = $c", legend = :outertopright, xlabel="β", ylabel="custo médio", size = (900, 500), 
            xtickfontsize = 20, ytickfontsize = 20, titlefontsize = 20, labelfontsize = 20, 
            linewidth = 3, legendfontsize=20, left_margin = 3mm, bottom_margin=6mm, right_margin = 2mm)
    plot!(β, [wait_and_return[1+k:4:end] for k in 0:4][index], label="Wait and Return", linewidth = 3)
    plot!(β, [naive_ignore[1+k:4:end] for k in 0:4][index], label="Naive Ignore", linewidth = 3)
    plot!(β, [naive_return[1+k:4:end] for k in 0:4][index], label="Naive Return", linewidth = 3)
    plot!(β, [compute_return[1+k:4:end] for k in 0:4][index], label="Compute Return", linewidth = 3)
end
get_plots_for_index(1)
for i in 1:4
    plt = get_plots_for_index(i)
    savefig(plt, "test/images/N=$N,capacity$(capacity[1:4][i]).pdf")
end

# Each algoritm per capacity
# Avg Cost vs. β for each capacity
function get_plots_for_index(index)
    algorithms = [wait_and_ignore, wait_and_return, naive_ignore, naive_return, compute_return]
    names = ["Wait and Ignore", "Wait and Return", "Naive Ignore", "Naive Return", "Compute Return"]
    β = [50, 100, 250, 500, 1000]
    algo_data = [algorithms[index][1+k:4:end] for k in 0:4]
    plot(β, algo_data[1], label="c=1", title="$(names[index])", 
        legend = :outertopright, xlabel="β", ylabel="custo médio", size = (600, 500), 
        xtickfontsize = 16, ytickfontsize = 16, titlefontsize = 20, labelfontsize = 20, 
        linewidth = 3, legendfontsize=20, left_margin = 5mm, bottom_margin=6mm)
    plot!(β, algo_data[2], label="c=2", linewidth = 3)
    plot!(β, algo_data[3], label="c=3", linewidth = 3)
    plot!(β, algo_data[4], label="c=5", linewidth = 3)
end
for i in 1:5
    plt = get_plots_for_index(i)
    savefig(plt, "test/images/N=$N,Algo $(["Wait and Ignore", "Wait and Return", "Naive Ignore", "Naive Return", "Compute Return"][i]).pdf")
end

# SMALL DATABASE RESULTS
data = XLSX.readdata("test/results_large.xlsx", "Sheet1", "A1:G171")
capacity = data[2:5:end, 2]
wait_and_ignore = Vector{Vector{Float64}}() #data[2:5:end, 5]
wait_and_return = Vector{Vector{Float64}}() #data[3:5:end, 5]
naive_ignore = Vector{Vector{Float64}}() #data[4:5:end, 5]
naive_return = Vector{Vector{Float64}}() #data[5:5:end, 5]
compute_return = Vector{Vector{Float64}}() #data[6:5:end, 5]

wait_and_ignore_i = Vector{Float64}()
wait_and_return_i = Vector{Float64}()
naive_ignore_i = Vector{Float64}()
naive_return_i = Vector{Float64}()
compute_return_i = Vector{Float64}()

i = 0
for c in capacity
    if c == 1
        push!(wait_and_ignore, wait_and_ignore_i)
        push!(wait_and_return, wait_and_return_i)
        push!(naive_ignore, naive_ignore_i)
        push!(naive_return, naive_return_i)
        push!(compute_return, compute_return_i)
        wait_and_ignore_i = []
        wait_and_return_i = []
        naive_ignore_i = []
        naive_return_i = []
        compute_return_i = []
    end
    
    push!(wait_and_ignore_i, data[2+i*5, 4])
    push!(wait_and_return_i, data[3+i*5, 4])
    push!(naive_ignore_i, data[4+i*5, 4])
    push!(naive_return_i, data[5+i*5, 4])
    push!(compute_return_i, data[6+i*5, 4])
    i += 1
end
push!(wait_and_ignore, wait_and_ignore_i)
push!(wait_and_return, wait_and_return_i)
push!(naive_ignore, naive_ignore_i)
push!(naive_return, naive_return_i)
push!(compute_return, compute_return_i)
wait_and_ignore = wait_and_ignore[2:end]
wait_and_return = wait_and_return[2:end]
naive_ignore = naive_ignore[2:end]
naive_return = naive_return[2:end]
compute_return = compute_return[2:end]

algo_cost(algo, index) = [row[index] for row in algo]

# Avg Cost vs. Capacity for each instance
function get_plots_for_index(index) 
    algo_names = ["Wait and Ignore", "Wait and Return", "Naive Ignore", "Naive Return", "Compute Return"]
    wait_and_ignore_cost = algo_cost(wait_and_ignore, index)
    wait_and_return_cost = algo_cost(wait_and_return, index)
    naive_ignore_cost = algo_cost(naive_ignore, index)
    naive_return_cost = algo_cost(naive_return, index)
    compute_return_cost = algo_cost(compute_return, index)
    costs = [wait_and_ignore_cost, wait_and_return_cost, naive_ignore_cost, naive_return_cost, compute_return_cost]
    min_costs = [minimum([costs[i][k] for i in 1:5]) for k in 1:length(wait_and_ignore_cost)]

    wait_and_ignore_wins = sum(wait_and_ignore_cost[i] ≈ min_costs[i] for i in 1:length(wait_and_ignore_cost))
    wait_and_return_wins = sum(wait_and_return_cost[i] ≈ min_costs[i] for i in 1:length(wait_and_ignore_cost))
    naive_ignore_wins = sum(naive_ignore_cost[i] ≈ min_costs[i] for i in 1:length(wait_and_ignore_cost))
    naive_return_wins = sum(naive_return_cost[i] ≈ min_costs[i] for i in 1:length(wait_and_ignore_cost))
    compute_return_wins = sum(compute_return_cost[i] ≈ min_costs[i] for i in 1:length(wait_and_ignore_cost))
    wins = [wait_and_ignore_wins, wait_and_return_wins, naive_ignore_wins, naive_return_wins, compute_return_wins]
    plt = bar(algo_names, wins, title="Capacidade = $(capacity[index])", legend=false, ylim=(0, 27),
                size = (900, 500), xtickfontsize = 14, ytickfontsize = 20, titlefontsize = 20, 
                labelfontsize = 20, legendfontsize=20, xrotation=7)#, left_margin = 5mm, bottom_margin=6mm)
    return plt
end
get_plots_for_index(4)
for i in 1:4
    plt = get_plots_for_index(i)
    savefig(plt, "test/images/wins_per_algo_c=$(capacity[i]).pdf")
end

function get_plots_for_index() 
    wait_and_ignore_cost = [mean(algo_cost(wait_and_ignore, index)) for index in 1:4]
    wait_and_return_cost = [mean(algo_cost(wait_and_return, index)) for index in 1:4]
    naive_ignore_cost = [mean(algo_cost(naive_ignore, index)) for index in 1:4]
    naive_return_cost = [mean(algo_cost(naive_return, index)) for index in 1:4]
    compute_return_cost = [mean(algo_cost(compute_return, index)) for index in 1:4]
    
    plt = plot(capacity[1:4], wait_and_ignore_cost, label="Wait and Ignore", 
                xlabel="capacidade", ylabel="custo médio", size = (900, 500), 
                xtickfontsize = 20, ytickfontsize = 20, titlefontsize = 20, labelfontsize = 20, 
                linewidth = 3, legendfontsize=20, left_margin = 5mm, bottom_margin=6mm)
    plot!(capacity[1:4], wait_and_return_cost, linewidth = 3, label="Wait and Return")
    plot!(capacity[1:4], naive_ignore_cost, linewidth = 3, label="Naive Ignore")
    plot!(capacity[1:4], naive_return_cost, linewidth = 3, label="Naive Return")
    plot!(capacity[1:4], compute_return_cost, linewidth = 3, label="Compute Return")
    return plt
end
plt = get_plots_for_index()
savefig(plt, "test/images/costs_db_small.pdf")

function get_plots_for_index(index)
    algo_names = ["Wait and Ignore", "Wait and Return", "Naive Ignore", "Naive Return", "Compute Return"]
    algos = [wait_and_ignore, wait_and_return, naive_ignore, naive_return, compute_return]
    algo_costs = vcat([algo_cost(algos[index], i) for i in 1:4]...)
    plt = violin(string.(capacity[1:4]), algo_costs, title=algo_names[index], legend=false,
                xlabel="capacidade", ylabel="custo médio", size = (600, 500), 
                xtickfontsize = 20, ytickfontsize = 16, titlefontsize = 20, labelfontsize = 20, 
                linewidth = 0, legendfontsize=20, left_margin = 2mm, bottom_margin=5mm)
    boxplot!(string.(capacity[1:4]), algo_costs, fillalpha=0.75, linewidth=2)
    return plt
end
get_plots_for_index(4)
for i in 1:5
    plt = get_plots_for_index(i)
    savefig(plt, "test/images/boxplot_db_small_$i.pdf")
end


# LARGE DB RESULTS
data = XLSX.readdata("test/results_large.xlsx", "Sheet1", "A1:G171")
capacity = data[2:5:end, 2]
wait_and_ignore = Vector{Vector{Any}}() #data[2:5:end, 5]
wait_and_return = Vector{Vector{Any}}() #data[3:5:end, 5]
naive_ignore = Vector{Vector{Any}}() #data[4:5:end, 5]
naive_return = Vector{Vector{Any}}() #data[5:5:end, 5]
compute_return = Vector{Vector{Any}}() #data[6:5:end, 5]

wait_and_ignore_i = Vector{Any}(missing, 5)
wait_and_return_i = Vector{Any}(missing, 5)
naive_ignore_i = Vector{Any}(missing, 5)
naive_return_i = Vector{Any}(missing, 5)
compute_return_i = Vector{Any}(missing, 5)

i = 0
for c in capacity
    if !ismissing(data[2+i*5, 1])
        push!(wait_and_ignore, wait_and_ignore_i)
        push!(wait_and_return, wait_and_return_i)
        push!(naive_ignore, naive_ignore_i)
        push!(naive_return, naive_return_i)
        push!(compute_return, compute_return_i)
        wait_and_ignore_i = Vector{Any}(missing, 5)
        wait_and_return_i = Vector{Any}(missing, 5)
        naive_ignore_i = Vector{Any}(missing, 5)
        naive_return_i = Vector{Any}(missing, 5)
        compute_return_i = Vector{Any}(missing, 5)
    end
    
    wait_and_ignore_i[c] = data[2+i*5, 4]
    wait_and_return_i[c] = data[3+i*5, 4]
    naive_ignore_i[c] = data[4+i*5, 4]
    naive_return_i[c] = data[5+i*5, 4]
    compute_return_i[c] = data[6+i*5, 4]
    i += 1
end
wait_and_ignore_i[c] = data[2+i*5, 4]
wait_and_return_i[c] = data[3+i*5, 4]
naive_ignore_i[c] = data[4+i*5, 4]
naive_return_i[c] = data[5+i*5, 4]
compute_return_i[c] = data[6+i*5, 4]

wait_and_ignore = wait_and_ignore[2:end]
wait_and_return = wait_and_return[2:end]
naive_ignore = naive_ignore[2:end]
naive_return = naive_return[2:end]
compute_return = compute_return[2:end]

algo_cost(algo, index) = filter(x -> !ismissing(x), [row[index] for row in algo])

# Avg Cost vs. Capacity for each instance
function get_plots_for_index(capacity) 
    algo_names = ["Wait and Ignore", "Wait and Return", "Naive Ignore", "Naive Return", "Compute Return"]
    wait_and_ignore_cost = algo_cost(wait_and_ignore, capacity)
    wait_and_return_cost = algo_cost(wait_and_return, capacity)
    naive_ignore_cost = algo_cost(naive_ignore, capacity)
    naive_return_cost = algo_cost(naive_return, capacity)
    compute_return_cost = algo_cost(compute_return, capacity)
    costs = [wait_and_ignore_cost, wait_and_return_cost, naive_ignore_cost, naive_return_cost, compute_return_cost]
    min_costs = [minimum([costs[i][k] for i in 1:5]) for k in 1:length(wait_and_ignore_cost)]

    wait_and_ignore_wins = sum(wait_and_ignore_cost[i] ≈ min_costs[i] for i in 1:length(wait_and_ignore_cost))
    wait_and_return_wins = sum(wait_and_return_cost[i] ≈ min_costs[i] for i in 1:length(wait_and_ignore_cost))
    naive_ignore_wins = sum(naive_ignore_cost[i] ≈ min_costs[i] for i in 1:length(wait_and_ignore_cost))
    naive_return_wins = sum(naive_return_cost[i] ≈ min_costs[i] for i in 1:length(wait_and_ignore_cost))
    compute_return_wins = sum(compute_return_cost[i] ≈ min_costs[i] for i in 1:length(wait_and_ignore_cost))
    wins = [wait_and_ignore_wins, wait_and_return_wins, naive_ignore_wins, naive_return_wins, compute_return_wins]
    plt = bar(algo_names, wins, title="Capacidade = $capacity", legend=false, ylim=(0, 12),
                size = (900, 500), xtickfontsize = 14, ytickfontsize = 20, titlefontsize = 20, 
                labelfontsize = 20, legendfontsize=20, xrotation=7)#, left_margin = 5mm, bottom_margin=6mm)
    return plt
end
get_plots_for_index(5)
for i in [1, 2, 3, 5]
    plt = get_plots_for_index(i)
    savefig(plt, "test/images/wins_per_algo_c=$i.pdf")
end

function get_plots_for_index() 
    wait_and_ignore_cost = [mean(algo_cost(wait_and_ignore, index)) for index in [1, 2, 3, 5]]
    wait_and_return_cost = [mean(algo_cost(wait_and_return, index)) for index in [1, 2, 3, 5]]
    naive_ignore_cost = [mean(algo_cost(naive_ignore, index)) for index in [1, 2, 3, 5]]
    naive_return_cost = [mean(algo_cost(naive_return, index)) for index in [1, 2, 3, 5]]
    compute_return_cost = [mean(algo_cost(compute_return, index)) for index in [1, 2, 3, 5]]
    
    plt = plot([1, 2, 3, 5], wait_and_ignore_cost, label="Wait and Ignore", 
                xlabel="capacidade", ylabel="custo médio", size = (900, 500), 
                xtickfontsize = 20, ytickfontsize = 20, titlefontsize = 20, labelfontsize = 20, 
                linewidth = 3, legendfontsize=20, left_margin = 5mm, bottom_margin=6mm)
    plot!([1, 2, 3, 5], wait_and_return_cost, linewidth = 3, label="Wait and Return")
    plot!([1, 2, 3, 5], naive_ignore_cost, linewidth = 3, label="Naive Ignore")
    plot!([1, 2, 3, 5], naive_return_cost, linewidth = 3, label="Naive Return")
    plot!([1, 2, 3, 5], compute_return_cost, linewidth = 3, label="Compute Return")
    return plt
end
plt = get_plots_for_index()
savefig(plt, "test/images/costs_db_small.pdf")

function get_plots_for_index(index)
    algo_names = ["Wait and Ignore", "Wait and Return", "Naive Ignore", "Naive Return", "Compute Return"]
    algos = [wait_and_ignore, wait_and_return, naive_ignore, naive_return, compute_return]
    algo_costs = vcat([algo_cost(algos[index], i) for i in [1, 2, 3, 5]]...)
    plt = violin(string.([1, 2, 3, 5]), algo_costs, title=algo_names[index], legend=false,
                xlabel="capacidade", ylabel="custo médio", size = (600, 500), 
                xtickfontsize = 20, ytickfontsize = 16, titlefontsize = 20, labelfontsize = 20, 
                linewidth = 0, legendfontsize=20, left_margin = 2mm, bottom_margin=5mm)
    boxplot!(string.([1, 2, 3, 5]), algo_costs, fillalpha=0.75, linewidth=2)
    return plt
end
get_plots_for_index(5)
for i in 1:5
    plt = get_plots_for_index(i)
    savefig(plt, "test/images/boxplot_db_small_$i.pdf")
end

# COMPETITIVE RATIO RESULTS
data = XLSX.readdata("test/competitive_ratio_worst5.xlsx", "Sheet1", "A1:H97")
capacity = data[2:6:end, 3]
data[2:25:end, 6]
wait_and_ignore = data[2:6:end, 5]
wait_and_return = data[3:6:end, 5]
naive_ignore = data[4:6:end, 5]
naive_return = data[5:6:end, 5]
compute_return = data[6:6:end, 5]

# Algorithms per β
# Avg Cost vs. Capacity for each β
function get_plots_for_index(index)
    β = [50, 250, 500, 1000][index]
    c = capacity[1:4]
    plt = plot(c, [wait_and_ignore[1+4*k:4+4*k] for k in 0:3][index], label="Wait and Ignore", 
                title="β = $β", xlabel="capacidade", ylabel="competitive ratio", size = (900, 500), 
                xtickfontsize = 20, ytickfontsize = 20, titlefontsize = 20, labelfontsize = 20, legend=:outertopright,
                linewidth = 3, legendfontsize=20, left_margin = 5mm, bottom_margin=6mm)
    plot!(c, [wait_and_return[1+4*k:4+4*k] for k in 0:3][index], label="Wait and Return", linewidth = 3)
    plot!(c, [naive_ignore[1+4*k:4+4*k] for k in 0:3][index], label="Naive Ignore", linewidth = 3)
    plot!(c, [naive_return[1+4*k:4+4*k] for k in 0:3][index], label="Naive Return", linewidth = 3)
    plot!(c, [compute_return[1+4*k:4+4*k] for k in 0:3][index], label="Compute Return", linewidth = 3)
    return plt
end
get_plots_for_index(1)
for i in 1:4
    plt = get_plots_for_index(i)
    savefig(plt, "test/images/CompetitiveRatio_beta$([50, 250, 500, 1000][i]).pdf")
end

# Algorithms per capacity
# Avg Cost vs. β for each capacity
function get_plots_for_index(index)
    c = capacity[1:4][index]
    β = [50, 250, 500, 1000]
    plot(β, [wait_and_ignore[1+k:4:end] for k in 0:3][index], label="Wait and Ignore", 
            title="c = $c", legend = :outertopright, xlabel="β", ylabel="competitive ratio", size = (900, 500), 
            xtickfontsize = 20, ytickfontsize = 20, titlefontsize = 20, labelfontsize = 20, 
            linewidth = 3, legendfontsize=20, left_margin = 5mm, bottom_margin=6mm, right_margin = 2mm)
    plot!(β, [wait_and_return[1+k:4:end] for k in 0:3][index], label="Wait and Return", linewidth = 3)
    plot!(β, [naive_ignore[1+k:4:end] for k in 0:3][index], label="Naive Ignore", linewidth = 3)
    plot!(β, [naive_return[1+k:4:end] for k in 0:3][index], label="Naive Return", linewidth = 3)
    plot!(β, [compute_return[1+k:4:end] for k in 0:3][index], label="Compute Return", linewidth = 3)
end
get_plots_for_index(4)
for i in 1:4
    plt = get_plots_for_index(i)
    savefig(plt, "test/images/CompetitiveRatio_capacity$(capacity[1:4][i]).pdf")
end

# Each algoritm per capacity
# Avg Cost vs. β for each capacity
function get_plots_for_index(index)
    algorithms = [wait_and_ignore, wait_and_return, naive_ignore, naive_return, compute_return]
    names = ["Wait and Ignore", "Wait and Return", "Naive Ignore", "Naive Return", "Compute Return"]
    β = [50, 250, 500, 1000]
    algo_data = [algorithms[index][1+k:4:end] for k in 0:3]
    plot(β, algo_data[1], label="c=1", title="$(names[index])", 
        legend = :outertopright, xlabel="β", ylabel="competitive ratio", size = (600, 500), 
        xtickfontsize = 16, ytickfontsize = 16, titlefontsize = 20, labelfontsize = 20, 
        linewidth = 3, legendfontsize=20, left_margin = 5mm, bottom_margin=6mm)
    plot!(β, algo_data[2], label="c=2", linewidth = 3)
    plot!(β, algo_data[3], label="c=3", linewidth = 3)
    plot!(β, algo_data[4], label="c=5", linewidth = 3)
end
get_plots_for_index(1)
for i in 1:5
    plt = get_plots_for_index(i)
    savefig(plt, "test/images/CompetitiveRatio_Algo $(["Wait and Ignore", "Wait and Return", "Naive Ignore", "Naive Return", "Compute Return"][i]).pdf")
end
