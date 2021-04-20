module WaitAndIgnore

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
        
        # TODO: Verify optimality - gets first elements of optimal route
        min_cost, end_t, opt_route = OfflineAlgorithm.run(graph, active_requests, capacity, cur_t)
        flat_route = collect(Iterators.flatten(opt_route))
        route_length = length(flat_route)
        travel_route = flat_route[1:min(route_length, capacity)]

        # TODO: Allow different requests with same destination
        find_req(dest) = findfirst(request -> request.destination == dest, active_requests)
        travel_requests = map(dest -> active_requests[find_req(dest)], travel_route)

        min_cost, end_t, opt_route = OfflineAlgorithm.run(graph, travel_requests, capacity, cur_t)
        cost += min_cost
        cur_t = end_t

        push!(route, travel_route)

        filter!(request -> !(request.destination in travel_route), remaining_requests)
    end

    return cost, cur_t, route
end

end # module