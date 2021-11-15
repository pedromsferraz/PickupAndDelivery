module GraphPreprocessing

using Graphs, SimpleWeightedGraphs, Combinatorics

# Used for preprocessing vertices optimal paths
# Maps a set of vertices to the path of minimum cost on the graph
opt_paths = Dict{Vector{Int64}, Vector{Int64}}()

# Distances between vertices
dists = Array{Float64, 2}(undef, 0, 2)

# Returns cost of path
# Assumes distances array has already been filled
function path_cost(path::Vector{Int}, initial_time::Float64=0.0)
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

function midway_path_cost(path::Vector{Int}, initial_time::Float64=0.0)
    cost = 0.0
    cur_t = initial_time
    current = path[1]
    for vertex in path[2:end]
        cur_t += dists[current, vertex]
        cost += cur_t
        current = vertex
    end
    cur_t += dists[current, 1]
    return cost, cur_t
end

# Finds the path that minimizes sum of times to reach each vertex
function best_path(vertices::Vector{Int}, 
                    initial_time::Float64=0.0)
    min_cost = Inf64
    end_time = Inf64
    opt_path = Vector{Int64}()

    s_vertices = sort(vertices)
    global opt_paths
    if haskey(opt_paths, s_vertices)
        min_cost, end_time = path_cost(opt_paths[s_vertices], initial_time)
        return min_cost, end_time, opt_paths[s_vertices]
    end

    for permutation in permutations(vertices)
        cost, time = path_cost(permutation, initial_time)
        if (cost < min_cost)
            min_cost = cost
            end_time = time # time is the same for every path
            opt_path = permutation
        end
    end

    opt_paths[s_vertices] = opt_path
    return min_cost, end_time, opt_path
end

# Finds the path that minimizes sum of times to reach each vertex
function greedy_best_path(vertices::Vector{Int}, initial_time::Float64=0.0)
    min_cost = 0.0
    cur_t = initial_time
    opt_path = Vector{Int64}()

    cur_vertex = 1
    while length(vertices) > 0
        cost = Inf64
        next_vertex = nothing
        for vertex in vertices
            if dists[cur_vertex, vertex] < cost
                cost = dists[cur_vertex, vertex]
                next_vertex = vertex
            end 
        end
        cur_t += cost
        min_cost += cur_t
        push!(opt_path, next_vertex)

        vertices = filter(vertex -> vertex != next_vertex, vertices)
        cur_vertex = next_vertex
    end

    s_vertices = sort(vertices)
    opt_paths[s_vertices] = opt_path
    return min_cost, cur_t, opt_path
end

function distance_graph(graph::SimpleWeightedGraph)
    g = SimpleWeightedGraph(nv(graph))
    dists = floyd_warshall_shortest_paths(graph, graph.weights).dists
    for i in 1:nv(graph)
        for j in i+1:nv(graph)
            add_edge!(g, i, j, dists[i, j])
        end
    end
    return g
end

function preprocess_dists(graph::SimpleWeightedGraph)
    global dists = floyd_warshall_shortest_paths(graph, graph.weights).dists
end

# Preprocesses graph with best path for every subset of vertices
# Subsequent queries only have to calculate the cost of the precomputed paths
# TODO: precomputed path cost query in O(1)
function preprocess(graph::SimpleWeightedGraph)
    global opt_paths = Dict()
    preprocess_dists(graph)
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
            cost, time = path_cost(permutation)
            if (cost < min_cost)
                min_cost = cost
                end_time = time # time is the same for every path
                opt_path = permutation
            end
        end

        opt_paths[subset] = opt_path
    end
end

end # module