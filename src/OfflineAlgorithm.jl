module OfflineAlgorithm

using ..DataModel, ..GraphPreprocessing, SimpleWeightedGraphs, Combinatorics

# Brute-force offline algorithm for a fixed requests permutation
function offline_algorithm(graph::SimpleWeightedGraph, 
                        requests::Vector{Request}, 
                        capacity::Int64, 
                        initial_t::Float64)
    N = length(requests)
    min_cost = Inf64
    end_t = initial_t
    opt_route = Vector{Vector{Int64}}()

    # Vector of possible release times
    release_times = map(request -> request.release_time, requests)
    push!(release_times, initial_t)
    release_times = unique(release_times)
    release_times = filter(time -> time >= initial_t, release_times)

    for release_time in release_times
        # Valid requests for a given release time
        valid_requests = filter(request -> request.release_time <= release_time, requests)

        # Get valid requests vertices
        N_valid = length(valid_requests)
        vertices = map(request -> request.destination, valid_requests)

        # For each possible size of group
        for k in 1:min(capacity, N_valid)
            cost = 0.0
            time = release_time
            route = Vector{Vector{Int64}}()

            # Wait and fulfill first k requests
            path = vertices[1:k]
            path_cost, time = GraphPreprocessing.path_cost(path, time)
            cost += path_cost
            push!(route, path)

            # Fulfill remaining requests recursively
            if (k+1 <= N)
                @inbounds remaining_requests = filter(request -> !(request in valid_requests[1:k]), requests)
                rec_cost, time, rec_route = offline_algorithm(graph, remaining_requests, capacity, time)
                cost += rec_cost
                route = vcat(route, rec_route)
            end

            # Update min cost of instance
            if (cost < min_cost)
                min_cost = cost
                end_t = time
                opt_route = route
            end
        end
    end

    return min_cost, end_t, opt_route
end

# Run brute-force optimal offline algorithm
function run(graph::SimpleWeightedGraph, 
                requests::Vector{Request}, 
                capacity::Int64,
                initial_t::Float64=0.0)
    min_cost = Inf64
    end_t = initial_t
    opt_route = Vector{Vector{Int64}}()
    
    GraphPreprocessing.preprocess_dists(graph)

    # Compute the best solution for each fixed permutation of requests
    for requests_perm in collect(permutations(requests))
        cost_perm, end_t_perm, route_perm = offline_algorithm(graph, requests_perm, capacity, initial_t)

        # Update the min cost between permutations
        if (cost_perm < min_cost)
            min_cost = cost_perm
            end_t = end_t_perm
            opt_route = route_perm
        end
    end


    return min_cost, end_t, opt_route
end

# Run brute-force optimal offline algorithm using Multithreading
function run_multithreaded(graph::SimpleWeightedGraph, 
                requests::Vector{Request}, 
                capacity::Int64,
                initial_t::Float64=0.0)
    min_cost_T = repeat([Inf64], Threads.nthreads())
    end_t_T = repeat([initial_t], Threads.nthreads())
    opt_route_T = Vector{Vector{Vector{Int64}}}(undef, Threads.nthreads())

    GraphPreprocessing.preprocess_dists(graph)

    # Compute the best solution for each fixed permutation of requests
    Threads.@threads for requests_perm in collect(permutations(requests))
        cost_perm, end_t_perm, route_perm = offline_algorithm(graph, requests_perm, capacity, initial_t)

        # Update the min cost between permutations in each thread
        if (cost_perm < min_cost_T[Threads.threadid()])
            min_cost_T[Threads.threadid()] = cost_perm
            end_t_T[Threads.threadid()] = end_t_perm
            opt_route_T[Threads.threadid()] = route_perm
        end
    end

    # Get the best solution among threads
    best = argmin(min_cost_T)
    @inbounds return min_cost_T[best], end_t_T[best], opt_route_T[best]
end

end # module