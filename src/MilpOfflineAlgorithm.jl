module MilpOfflineAlgorithm

# Gurobi setup
# ENV["GUROBI_HOME"] = "/Library/gurobi912/mac64"
# using Pkg
# Pkg.add("Gurobi")
# Pkg.build("Gurobi")
using ..DataModel, ..GraphPreprocessing, LightGraphs, SimpleWeightedGraphs, JuMP, Gurobi

# Mixed-integer programming implementation of offline algorithm
# Note: This version assume all requests are already active
function offline_algorithm(graph::SimpleWeightedGraph, 
                            requests::Vector{Request}, 
                            capacity::Int64, 
                            initial_t::Float64)
    g = copy(graph)
    Kmax = ceil(Int, length(requests) / capacity)
    N = nv(g)

    # add n+1 vertex
    add_vertex!(g)
    for adj in neighbors(g, 1)
        add_edge!(g, adj, N+1, g.weights[1, adj])
    end
    GraphPreprocessing.preprocess_dists(g)

    # complete g with distances
    for i in 1:N+1
        for j in 1:N+1
            if !has_edge(g, i, j) && i != j
                add_edge!(g, i, j, GraphPreprocessing.dists[i, j])
            end
        end
    end
    add_edge!(g, 1, N+1, 1e-100)
    GraphPreprocessing.preprocess_dists(g)

    # d[i] = 1 if there is demand on node i, 0 otherwise
    d = zeros(N+1)
    for request in requests
        d[request.destination] = 1
    end
    
    model = Model(Gurobi.Optimizer)

    # x[i, j, k] -> if edge (i, j) is being used in route k
    @variable(model, x[1:N+1, 1:N+1, 1:Kmax], Bin)

    # ω[i, k] -> time at which request node i is served in route k
    @variable(model, ω[1:N+1, 1:Kmax] >= 0.0)

    # objective: minimize sum of serving time for every request
    @objective(model, Min, sum(d[i] * sum(ω[i, k] for i in 2:N) for k in 1:Kmax))

    # each request i must be assigned to a unique route k
    for i in 2:N
        adj = outneighbors(g, i)
        if d[i] == 1
            @constraint(model, sum(sum(x[i, j, k] for j in adj) for k in 1:Kmax) == 1)
        end
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
        adj = inneighbors(g, 1)
        @constraint(model, sum(x[i, N+1, k] for i in adj) == 1)
    end

    # schedule constraints
    for k in 1:Kmax
        for edge in edges(g)
            i = src(edge)
            j = dst(edge)
            # t_ij = weight(edge)
            t_ij = GraphPreprocessing.dists[i, j]
            @constraint(model, x[i, j, k] => { ω[i, k] + t_ij - ω[j, k] <= 0 })
            @constraint(model, x[j, i, k] => { ω[j, k] + t_ij - ω[i, k] <= 0 })
        end
    end

    # capacity constraints
    for k in 1:Kmax
        @constraint(model, sum(d[i] * sum(x[i, j, k] for j in outneighbors(g, i)) for i in 2:N) <= capacity)
    end

    # single vehicle constraint
    for k in 1:Kmax-1
        @constraint(model, ω[1, k+1] == ω[N+1, k])
    end

    # initial time constraint
    @constraint(model, ω[1, 1] == initial_t)

    optimize!(model)

    termination_status(model)
    min_cost = objective_value(model)
    end_t = value(ω[N+1, Kmax])
    value.(x)
    value.(ω)

    return min_cost, end_t
end

function run(graph::SimpleWeightedGraph, 
                requests::Vector{Request}, 
                capacity::Int64,
                initial_t::Float64=0.0)
    return offline_algorithm(graph, requests, capacity, initial_t)
end

end # module