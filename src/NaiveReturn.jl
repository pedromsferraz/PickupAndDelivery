module NaiveReturn

using ..DataModel, ..GraphPreprocessing, ..OfflineAlgorithm, ..MilpOfflineAlgorithm, SimpleWeightedGraphs

const USES_MILP = true

function call_offline_algorithm(graph, requests, capacity, cur_t, request_limit)
    if request_limit == 0 || length(requests) <= request_limit
        return OfflineAlgorithm.run_multithreaded(graph, requests, capacity, cur_t)
    else
        return MilpOfflineAlgorithm.run(graph, requests, capacity, cur_t, full_search=false, time_limit=60.0)
    end
end

function run(graph::SimpleWeightedGraph, 
            requests::Vector{Request}, 
            capacity::Int64,
            request_limit::Int64 = 0)
    GraphPreprocessing.preprocess_dists(graph)
    
    N = length(requests)
    
    cost = 0.0
    cur_t = 0.0
    route = Vector{Vector{Int64}}()
    remaining_requests = copy(requests)

    while !isempty(remaining_requests)
        first_valid = minimum(map(request -> request.release_time, remaining_requests))
        cur_t = max(cur_t, first_valid)
        valid_requests = filter(request -> request.release_time <= cur_t, remaining_requests)
        
        if !USES_MILP && request_limit != 0 && length(valid_requests) > request_limit
            throw(RequestLimit("Wait and Return: limite de pedidos excedido - tentou executar o algoritmo offline com $(length(valid_requests)) pedidos ativos."))
        end

        # Calculate optimal route to serve valid requests
        min_cost, end_t, opt_route = call_offline_algorithm(graph, valid_requests, capacity, cur_t, request_limit)

        # Calculate requests that will become valid during the route and sort by valid time
        partway_valid_requests = filter(request -> cur_t < request.release_time < end_t, remaining_requests)
        sort!(partway_valid_requests, by = request -> request.release_time)

        # Calculate maximal distance from requests to origin
        requests_dests = map(request -> request.destination, valid_requests)
        maximimal_dist = maximum(vertex -> GraphPreprocessing.dists[1, vertex], requests_dests)

        flat_route = opt_route[1] # collect(Iterators.flatten(opt_route))
        route_length = length(flat_route)
        flat_route = flat_route[1:min(route_length, capacity)]
        travel_route = Vector{Int64}()
        prev_vertex = 1
        num_new_requests = 0
        num_unserved_requests = length(flat_route)
        should_return_origin = false
        for vertex in flat_route
            next_t = cur_t + GraphPreprocessing.dists[prev_vertex, vertex]

            # Check if request(s) will become valid during this edge
            partedge_valid_requests = filter(request -> cur_t < request.release_time < next_t, partway_valid_requests)
            sort!(partedge_valid_requests, by = request -> request.release_time)

            for partedge_request in partedge_valid_requests
                return_prev_time = partedge_request.release_time - cur_t
                return_origin_time = return_prev_time + GraphPreprocessing.dists[prev_vertex, 1]
                num_new_requests += 1
                
                # If condition is met, return to origin
                if return_origin_time / maximimal_dist <= num_new_requests / (num_new_requests + num_unserved_requests)
                    # Time travelled so far + time to return to vertex
                    cur_t += 2*return_prev_time
                    should_return_origin = true
                    push!(travel_route, 0)
                    break
                end
            end

            if should_return_origin
                break
            end
            cur_t += GraphPreprocessing.dists[prev_vertex, vertex]
            cost += cur_t
            prev_vertex = vertex
            num_unserved_requests -= 1
            push!(travel_route, vertex)
        end
        cur_t += GraphPreprocessing.dists[prev_vertex, 1]
                
        push!(route, travel_route)
        filter!(request -> !(request.destination in travel_route), remaining_requests)
    end

    return cost, cur_t, route
end

end # module