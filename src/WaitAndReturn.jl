module WaitAndReturn

using ..DataModel, ..GraphPreprocessing, ..OfflineAlgorithm, SimpleWeightedGraphs

function run(graph::SimpleWeightedGraph, 
            requests::Vector{Request}, 
            capacity::Int64)
    GraphPreprocessing.preprocess(graph)
    
    N = length(requests)
    active_time(request) = max(request.release_time, GraphPreprocessing.dists[1, request.destination])
    
    cost = 0.0
    cur_t = 0.0
    route = Vector{Vector{Int64}}()
    remaining_requests = copy(requests)

    while !isempty(remaining_requests)
        first_active = minimum(map(active_time, remaining_requests))
        cur_t = max(cur_t, first_active)
        active_requests = filter(request -> active_time(request) <= cur_t, remaining_requests)
        
        if length(active_requests) >= capacity
            # Serve first elements of optimal route
            min_cost, end_t, opt_route = OfflineAlgorithm.run(graph, active_requests, capacity, cur_t)
            flat_route = collect(Iterators.flatten(opt_route))
            route_length = length(flat_route)
            travel_route = flat_route[1:capacity]

            # TODO: Allow different requests with same destination
            find_req(dest) = findfirst(request -> request.destination == dest, active_requests)
            travel_requests = map(dest -> active_requests[find_req(dest)], travel_route)

            min_cost, end_t, opt_route = OfflineAlgorithm.run(graph, travel_requests, capacity, cur_t)
            cost += min_cost
            cur_t = end_t

            push!(route, travel_route)
            filter!(request -> !(request.destination in travel_route), remaining_requests)
        else
            # Calculate optimal route to serve active requests
            min_cost, end_t, opt_route = OfflineAlgorithm.run(graph, active_requests, capacity, cur_t)

            # Calculate requests that will become active during the route and sort by active time
            partway_active_requests = filter(request -> cur_t < active_time(request) < end_t, remaining_requests)
            sort!(partway_active_requests, by = request -> active_time(request))

            # Calculate maximal distance from requests to origin
            requests_dests = map(request -> request.destination, active_requests)
            maximimal_dist = maximum(vertex -> GraphPreprocessing.dists[1, vertex], requests_dests)

            flat_route = collect(Iterators.flatten(opt_route))
            travel_route = Vector{Int64}()
            prev_vertex = 1
            num_new_requests = 0
            num_unserved_requests = length(active_requests)
            should_return_origin = false
            for vertex in flat_route
                next_t = cur_t + GraphPreprocessing.dists[prev_vertex, vertex]

                # Check if request(s) will become active during this edge
                partedge_active_requests = filter(request -> cur_t < active_time(request) < next_t, partway_active_requests)
                sort!(partedge_active_requests, by = request -> active_time(request))

                for partedge_request in partedge_active_requests
                    return_prev_time = active_time(partedge_request) - cur_t
                    return_origin_time = return_prev_time + GraphPreprocessing.dists[prev_vertex, 1]
                    num_new_requests += 1
                    
                    # If condition is met, return to origin
                    if return_origin_time / maximimal_dist <= num_new_requests / (num_new_requests + num_unserved_requests)
                        # Time travelled so far + time to return to vertex
                        cur_t += 2*return_prev_time
                        should_return_origin = true
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
    end

    return cost, cur_t, route
end

end # module