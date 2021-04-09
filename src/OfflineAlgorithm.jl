module OfflineAlgorithm

using ..DataModel, LightGraphs, SimpleWeightedGraphs, Combinatorics

# Used for preprocessing vertices optimal paths
# Maps a set of vertices to the path of minimum cost on the graph
opt_paths = Dict{Vector{Int64}, Vector{Int64}}()

# Distances between vertices
dists = Array{Float64, 2}(undef, 0, 2)

# Last graph ran by the algorithm
global_graph = SimpleWeightedGraph()

# Returns cost of path
function path_cost(graph::SimpleWeightedGraph, 
                    path::Vector{Int}, 
                    initial_time::Float64=0.0)
    cost = 0.0
    cur_t = initial_time
    current = 1
    for vertex in path
        cur_t += dists[current, vertex]
        cost += cur_t
        current = vertex
    end
    cur_t += dists[current, 1]
    return cost, cur_t
end

# Finds the path that minimizes sum of times to reach each vertex
function best_path(graph::SimpleWeightedGraph, 
                    vertices::Vector{Int}, 
                    initial_time::Float64=0.0)
    min_cost = Inf64
    end_time = Inf64
    opt_path = Vector{Int64}()

    s_vertices = sort(vertices)
    global opt_paths
    if haskey(opt_paths, s_vertices)
        min_cost, end_time = path_cost(graph, opt_paths[s_vertices], initial_time)
        return min_cost, end_time, opt_paths[s_vertices]
    end

    for permutation in permutations(vertices)
        cost, time = path_cost(graph, permutation, initial_time)
        if (cost < min_cost)
            min_cost = cost
            end_time = time # time is the same for every path
            opt_path = permutation
        end
    end

    opt_paths[s_vertices] = opt_path
    return min_cost, end_time, opt_path
end

# Preprocesses graph with best path for every subset of vertices
# Subsequent queries only have to calculate the cost of the precomputed paths
# TODO: precomputed path cost query in O(1)
function precompute_best_paths(graph::SimpleWeightedGraph)
    global opt_paths = Dict()
    global dists = floyd_warshall_shortest_paths(graph, graph.weights).dists
    N = nv(graph)
    vertices = [2:N...]
    
    for subset in powerset(vertices)
        if length(subset) < 2
            continue
        end

        min_cost = Inf64
        end_time = Inf64
        opt_path = Vector{Int64}()
        
        for permutation in permutations(subset)
            cost, time = path_cost(graph, permutation)
            if (cost < min_cost)
                min_cost = cost
                end_time = time # time is the same for every path
                opt_path = permutation
            end
        end

        opt_paths[subset] = opt_path
    end

    global global_graph = graph
end

# Brute-force offline optimal algorithm recursive function
function offline_algorithm(graph::SimpleWeightedGraph, 
                        requests::Vector{Request}, 
                        capacity::Int64, 
                        initial_t::Float64=0.0)
    N = length(requests)
    min_cost = Inf64
    end_t = initial_t
    opt_route = Vector{Vector{Int}}()

    # Vector of possible release times
    release_times = map(request -> request.release_time, requests)
    push!(release_times, initial_t)
    release_times = unique(release_times)
    release_times = filter(time -> time >= initial_t, release_times)

    for release_time in release_times
        # Valid requests for a given release time
        valid_requests = filter(request -> request.release_time <= release_time, requests)

        # Get best permutation of valid requests
        for requests_perm in permutations(valid_requests)
            cost_perm = Inf64
            end_t_perm = release_time
            route_perm = Vector{Vector{Int}}()
            N_perm = length(requests_perm)

            vertices = map(request -> request.destination, requests_perm)

            # For each possible size of group
            for k in 1:min(capacity, N_perm)
                cost = 0.0
                time = release_time
                route = Vector{Vector{Int}}()

                # Wait and fulfill first k requests
                ful_vert = vertices[1:k]
                path_cost, time, path = best_path(graph, ful_vert, time)
                cost += path_cost
                push!(route, path)

                # Fulfill remaining requests recursively
                if (k+1 <= N)
                    remaining_requests = filter(request -> !(request in requests_perm[1:k]), requests)
                    rec_cost, time, rec_route = offline_algorithm(graph, remaining_requests, capacity, time)
                    cost += rec_cost
                    route = vcat(route, rec_route)
                end

                # Update min cost of permutation
                if (cost < cost_perm)
                    cost_perm = cost
                    end_t_perm = time
                    route_perm = route
                end
            end

            # Update min cost of instance
            if (cost_perm < min_cost)
                min_cost = cost_perm
                end_t = end_t_perm
                opt_route = route_perm
            end
        end
    end

    return min_cost, end_t, opt_route
end

# Run brute-force optimal offline algorithm
function run(graph::SimpleWeightedGraph, 
                requests::Vector{Request}, 
                capacity::Int64)
    global global_graph
    if graph != global_graph
        precompute_best_paths(graph)
    end
    return offline_algorithm(graph, requests, capacity)
end

end # module