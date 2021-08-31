module ComputeReturn

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

        flat_route = opt_route[1] # collect(Iterators.flatten(opt_route))
        route_length = length(flat_route)
        flat_route = flat_route[1:min(route_length, capacity)]
        travel_route = Vector{Int64}()
        prev_vertex = 1
        should_return_origin = false
        for (i, vertex) in enumerate(flat_route)
            next_t = cur_t + GraphPreprocessing.dists[prev_vertex, vertex]

            # Check if request(s) will become valid during this edge
            partedge_valid_requests = filter(request -> cur_t < request.release_time < next_t, partway_valid_requests)
            sort!(partedge_valid_requests, by = request -> request.release_time)

            for partedge_request in partedge_valid_requests
                # Time at which the decision to go back to the origin is being made
                decision_time = partedge_request.release_time
                return_prev_time = decision_time - cur_t                
                released_remaining_requests = filter(request -> request.release_time <= decision_time, remaining_requests)
                waiting_requests = filter(request -> !(request.destination in flat_route), released_remaining_requests)

                # Cost of continuing in this route
                if i == 1
                    current_vertex_index = 1
                    remaining_route = [1, flat_route[1:end]...]
                else
                    current_vertex_index = i-1
                    remaining_route = flat_route[current_vertex_index:end]
                end
                continue_route_initial_cost, continue_route_end_t = GraphPreprocessing.midway_path_cost(remaining_route, cur_t)
                if length(waiting_requests) > 0
                    continue_route_remaning_cost, _, _ = call_offline_algorithm(graph, waiting_requests, capacity, continue_route_end_t, request_limit)
                else 
                    continue_route_remaning_cost = 0
                end
                continue_route_cost = continue_route_initial_cost + continue_route_remaning_cost

                # Cost of returning to origin
                origin_arrival_time = decision_time + return_prev_time + GraphPreprocessing.dists[prev_vertex, 1]
                if i == 1
                    origin_remaining_requests = released_remaining_requests
                else
                    origin_remaining_requests = filter(request -> !(request.destination in flat_route[1:current_vertex_index]), released_remaining_requests)
                end
                return_route_cost, _, _ = call_offline_algorithm(graph, origin_remaining_requests, capacity, origin_arrival_time, request_limit)
                
                # If condition is met, return to origin
                if return_route_cost < continue_route_cost
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
            push!(travel_route, vertex)
        end
        cur_t += GraphPreprocessing.dists[prev_vertex, 1]
    
        push!(route, travel_route)
        filter!(request -> !(request.destination in travel_route), remaining_requests)
    end

    return cost, cur_t, route
end

end # module