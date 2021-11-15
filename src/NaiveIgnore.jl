module NaiveIgnore

using ..DataModel, ..GraphPreprocessing, ..OfflineAlgorithm, ..MilpOfflineAlgorithm, SimpleWeightedGraphs

const USES_MILP = true

function call_offline_algorithm(graph, requests, capacity, cur_t, request_limit)
    if request_limit == 0 || length(requests) <= request_limit
        return OfflineAlgorithm.run_multithreaded(graph, requests, capacity, cur_t)
    else
        return MilpOfflineAlgorithm.run(graph, requests, capacity, cur_t, full_search=false, time_limit=20.0)
    end
end

function run(graph::SimpleWeightedGraph, 
            requests::Vector{Request}, 
            capacity::Int64,
            request_limit::Int64 = 0,
            milp_limit::Int64 = 0)
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
        
        if (!USES_MILP && request_limit != 0 && length(valid_requests) > request_limit) || (USES_MILP && milp_limit != 0 && length(valid_requests) > milp_limit)
            throw(RequestLimit("Naive Ignore: limite de pedidos excedido - tentou executar o algoritmo offline com $(length(valid_requests)) pedidos ativos."))
        end

        min_cost, end_t, opt_route = call_offline_algorithm(graph, valid_requests, capacity, cur_t, request_limit)
        flat_route = opt_route[1] #collect(Iterators.flatten(opt_route))
        route_length = length(flat_route)
        travel_route = flat_route[1:min(route_length, capacity)]

        # TODO: Allow different requests with same destination
        find_req(dest) = findfirst(request -> request.destination == dest, valid_requests)
        travel_requests = map(dest -> valid_requests[find_req(dest)], travel_route)

        min_cost, end_t, opt_route = call_offline_algorithm(graph, travel_requests, capacity, cur_t, request_limit)
        cost += min_cost
        cur_t = end_t

        push!(route, travel_route)

        filter!(request -> !(request.destination in travel_route), remaining_requests)
    end

    return cost, cur_t, route
end

end # module