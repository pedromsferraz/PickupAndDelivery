module MilpOfflineAlgorithm

# Gurobi setup
# ENV["GUROBI_HOME"] = "/Library/gurobi912/mac64/"
# using Pkg
# Pkg.add("Gurobi")
# Pkg.build("Gurobi")
using ..DataModel, ..GraphPreprocessing, Graphs, SimpleWeightedGraphs, JuMP, Gurobi

function get_route(x::Array{Float64,3}, d::Vector{Int64})
    subgraph_route = Vector{Vector{Int64}}()
    N = size(x, 1)
    Kmax = size(x, 3)

    original_vertices = findall(x -> x == 1, d)
    for k in 1:Kmax
        k_route = Vector{Int64}()
        i = 1
        while i != N
            i = findfirst(e -> e ≈ 1.0, value.(x[i, :, k]))
            if i != N
                push!(k_route, i)
            end
        end
        if length(k_route) > 0
            push!(subgraph_route, k_route)
        end
    end

    recover_vertex(i::Int) = original_vertices[i]
    route = map(i -> recover_vertex.(i .- 1), subgraph_route)

    return route, subgraph_route
end

function get_end_t(g, subgraph_route, ω, Kmax)
    subgraph_last_vertex = collect(Iterators.flatten(subgraph_route))[end]
    return value(ω[subgraph_last_vertex, Kmax]) + g.weights[subgraph_last_vertex, 1]
end

function preprocess_graph(graph::SimpleWeightedGraph, d::Vector{Int64})
    g = copy(graph)
    N = nv(graph)

    # add n+1 vertex
    add_vertex!(g)
    for adj in neighbors(g, 1)
        add_edge!(g, adj, N+1, g.weights[1, adj])
    end
    add_edge!(g, 1, N+1, 1e-100)
    g = GraphPreprocessing.distance_graph(g)

    # remove unused vertices
    for i in Iterators.reverse(2:N)
        if d[i] == 0
            rem_vertex!(g, i)
        end
    end
    
    return g
end

function get_next_vertex(g, current_vertex, request_dest)
    w = filter(x -> x != 0.0, collect(g.weights[current_vertex, request_dest]))
    minw = minimum(w)
    argminw = argmin(w)
    return minw, request_dest[argminw]
end

function greedy_warm_start(graph::SimpleWeightedGraph, 
                        requests::Vector{Request}, 
                        capacity::Int64, 
                        initial_t::Float64)
    Kmin = ceil(Int, length(requests) / capacity)
    N = nv(graph)

    # d[i] = 1 if there is demand on node i, 0 otherwise
    d = zeros(Int, N)
    for request in requests
        d[request.destination] = 1
    end
    g = preprocess_graph(graph, d)
    N = nv(g) - 1
    x = zeros(N+1, N+1, Kmin)
    ω = zeros(N+1, Kmin)

    request_dest = [i for i in 2:N]
    cur_vertex = 1
    cur_t = initial_t
    cur_k = 1

    ω[cur_vertex, cur_k] = cur_t
    while length(request_dest) > 0
        cur_vertex = 1
        fulfilled_requests = 0
        while length(request_dest) > 0 && fulfilled_requests < capacity
            cost, next_vertex = get_next_vertex(g, cur_vertex, request_dest)
            cur_t += cost
            x[cur_vertex, next_vertex, cur_k] = 1.0
            ω[next_vertex, cur_k] = cur_t

            request_dest = filter(vertex -> vertex != next_vertex, request_dest)
            cur_vertex = next_vertex
            fulfilled_requests += 1
        end
        next_vertex = N+1
        cur_t += g.weights[cur_vertex, N+1]
        x[cur_vertex, next_vertex, cur_k] = 1.0
        ω[next_vertex, cur_k] = cur_t
        cur_k += 1
        if cur_k <= Kmin
            ω[1, cur_k] = cur_t
        end
    end

    return x, ω
end

# Mixed-integer programming implementation of offline algorithm
# Note: This version assume all requests are already active
function offline_algorithm(graph::SimpleWeightedGraph, 
                            requests::Vector{Request}, 
                            capacity::Int64, 
                            initial_t::Float64,
                            full_search::Bool,
                            time_limit::Float64,
                            warm_start::Bool)
    Kmin = ceil(Int, length(requests) / capacity)
    N = nv(graph)

    # d[i] = 1 if there is demand on node i, 0 otherwise
    d = zeros(Int, N)
    for request in requests
        d[request.destination] = 1
    end
    g = preprocess_graph(graph, d)
    N = nv(g) - 1

    min_cost = Inf
    end_t = Inf
    opt_route = Vector{Vector{Int64}}()

    Klist = full_search ? (Kmin:length(requests)) : [Kmin]
    # iterate through all possible 
    for Kmax in Klist
        model = Model(Gurobi.Optimizer)
        set_optimizer_attribute(model, "TimeLimit", time_limit)
        set_optimizer_attribute(model, "NumericFocus", 2)
        # set_silent(model)

        # x[i, j, k] -> if edge (i, j) is being used in route k
        @variable(model, x[1:N+1, 1:N+1, 1:Kmax], Bin)

        # ω[i, k] -> time at which request node i is served in route k
        @variable(model, ω[1:N+1, 1:Kmax] >= 0.0)

        # set warm start
        if warm_start
            initial_x, initial_ω = greedy_warm_start(graph, requests, capacity, initial_t)
            set_start_value.(x, initial_x)
            set_start_value.(ω, initial_ω)
        end

        # objective: minimize sum of serving time for every request
        @objective(model, Min, sum(sum(ω[i, k] for k in 1:Kmax) for i in 2:N))

        # each request i must be assigned to a unique route k
        for i in 2:N
            adj = outneighbors(g, i)
            @constraint(model, sum(sum(x[i, j, k] for j in adj) for k in 1:Kmax) == 1)
        end

        # flow constraints
        for k in 1:Kmax
            adj = outneighbors(g, 1)
            @constraint(model, sum(x[1, j, k] for j in adj) == 1)
        end

        for k in 1:Kmax
            for j in 2:N
                adj_in = inneighbors(g, j)
                adj_out = outneighbors(g, j)
                @constraint(model, sum(x[i, j, k] for i in adj_in) - sum(x[j, i, k] for i in adj_out) == 0)
            end
        end

        for k in 1:Kmax
            adj = inneighbors(g, N+1)
            @constraint(model, sum(x[i, N+1, k] for i in adj) == 1)
        end

        # schedule constraints
        for k in 1:Kmax
            for edge in edges(g)
                i = src(edge)
                j = dst(edge)
                t_ij = weight(edge)
                @constraint(model, x[i, j, k] => { ω[i, k] + t_ij - ω[j, k] <= 0 })
                @constraint(model, x[j, i, k] => { ω[j, k] + t_ij - ω[i, k] <= 0 })
            end
        end

        # capacity constraints
        for k in 1:Kmax
            @constraint(model, sum(sum(x[i, j, k] for j in outneighbors(g, i)) for i in 2:N) <= capacity)
        end

        # single vehicle constraint
        for k in 1:Kmax-1
            @constraint(model, ω[1, k+1] == ω[N+1, k])
        end

        # initial time constraint
        @constraint(model, ω[1, 1] == initial_t)

        optimize!(model)
        
        if (time_limit < Inf64 || termination_status(model) == MOI.OPTIMAL) &&
                !isapprox(objective_value(model), min_cost, atol=1e-3) && 
                objective_value(model) < min_cost
            min_cost = objective_value(model)
            opt_route, subgraph_route = get_route(value.(x), d)
            end_t = get_end_t(g, subgraph_route, ω, Kmax)
        end
    end

    return min_cost, end_t, opt_route
end

function run(graph::SimpleWeightedGraph, 
                requests::Vector{Request}, 
                capacity::Int64,
                initial_t::Float64=0.0;
                full_search::Bool=true,
                time_limit::Float64=Inf64,
                warm_start::Bool=true)
    return offline_algorithm(graph, requests, capacity, initial_t, full_search, time_limit, warm_start)
end

end # module