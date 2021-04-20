module OfflineAlgorithm

using ..DataModel, ..Graph, SimpleWeightedGraphs, Combinatorics

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
                path_cost, time, path = Graph.best_path(ful_vert, time)
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
    Graph.preprocess(graph)
    return offline_algorithm(graph, requests, capacity)
end

end # module