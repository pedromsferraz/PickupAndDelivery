using Test, FinalProjectPedroFerraz
using LightGraphs, SimpleWeightedGraphs, StatsBase, Random #, Plots

@testset "Offline algorithm tests" begin
    @testset "Simple small graph test" begin
        # Create simple weighted graph
        srcs = [1, 1, 2]
        dsts = [2, 3, 3]
        wgts = [3., 10., 1.]
        g = SimpleWeightedGraph(srcs, dsts, wgts);

        # Define requests
        N_req = 2
        release_times = [5, 40]
        request_vertices = [2, 3]
        requests = map(i -> DataModel.Request(release_times[i], request_vertices[i]), 1:N_req)

        # Define capacity
        capacity = 2

        # Run offline algorithm tests
        min_cost, end_t, opt_route = OfflineAlgorithm.run(g, requests, capacity)
        @test min_cost ≈ 52.0
        @test end_t ≈ 48.0
        @test opt_route == [[2], [3]]

        requests[2].release_time = 10
        min_cost, end_t, opt_route = OfflineAlgorithm.run(g, requests, capacity)
        @test min_cost ≈ 23.0
        @test end_t ≈ 19.0
        @test opt_route == [[2], [3]]

        requests[1].release_time = 10
        requests[2].release_time = 50
        min_cost, end_t, opt_route = OfflineAlgorithm.run(g, requests, capacity)
        @test min_cost ≈ 67.0
        @test end_t ≈ 58.0
        @test opt_route == [[2], [3]]
    end

    @testset "Unordered requests" begin
        # Create simple weighted graph
        srcs = [1, 1, 2]
        dsts = [2, 3, 3]
        wgts = [1., 100., 10000.]
        g = SimpleWeightedGraph(srcs, dsts, wgts);

        # Define requests
        N_req = 2
        release_times = [1, 0]
        request_vertices = [2, 3]
        requests = map(i -> DataModel.Request(release_times[i], request_vertices[i]), 1:N_req)

        # Define capacity
        capacity = 2

        # Run offline algorithm tests
        min_cost, end_t, opt_route = OfflineAlgorithm.run(g, requests, capacity)
        @test min_cost ≈ 105.0
        @test end_t ≈ 203.0
        @test opt_route == [[2], [3]]
    end

    @testset "Big edge weight" begin
        # Create simple weighted graph
        srcs = [1, 1, 2]
        dsts = [2, 3, 3]
        wgts = [1., 1000., 1.]
        g = SimpleWeightedGraph(srcs, dsts, wgts);

        # Define requests
        N_req = 2
        release_times = [0, 1000]
        request_vertices = [2, 3]
        requests = map(i -> DataModel.Request(release_times[i], request_vertices[i]), 1:N_req)

        # Define capacity
        capacity = 2

        # Run offline algorithm tests
        min_cost, end_t, opt_route = OfflineAlgorithm.run(g, requests, capacity)
        @test min_cost ≈ 1003.0
        @test end_t ≈ 1004.0
        @test opt_route == [[2], [3]]
    end

    @testset "Jamboard example" begin
        # Define and generate N_points in R^2
        N_points = 10
        points = [0 1 2 3 4 0 10 20 30 40; 10 10 10 10 10 -50 -50 -50 -50 -50]
        
        # Create origin -> origin will be vertex 1
        points = hcat([0, 0], points)

        # Plot points for verification
        # plot = Plots.plot(points[1, :], points[2, :], seriestype = :scatter, legend = false)
        # mkpath("images")
        # Plots.savefig(plot, "images/points.pdf")

        # Create simple weighted graph from euclidean graph
        g, dists = euclidean_graph(points, p=2)
        srcs = Vector{Int}()
        dsts = Vector{Int}()
        wgts = Vector{Float64}()
        for e in dists
            push!(srcs, e.first.src)
            push!(dsts, e.first.dst)
            push!(wgts, e.second)
        end
        g = SimpleWeightedGraph(srcs, dsts, wgts);

        # Define requests with equal release time
        N_req = 5
        release_times = zeros(N_req)
        request_vertices = sample(2:nv(g), N_req ,replace=false)
        requests = map(i -> DataModel.Request(release_times[i], request_vertices[i]), 1:N_req)

        # Run offline algorithm tests
        capacity = 5
        min_cost, end_t, opt_route = OfflineAlgorithm.run(g, requests, capacity)
        @test opt_route == [sort(request_vertices)]

        capacity = 1
        min_cost, end_t, opt_route = OfflineAlgorithm.run(g, requests, capacity)
        @test opt_route == map(el -> [el], sort(request_vertices))

        # Define requests with out of order release times
        release_times = shuffle([1:N_req...] .* 0.0001)
        request_vertices = sample(2:nv(g), N_req ,replace=false)
        requests = map(i -> DataModel.Request(release_times[i], request_vertices[i]), 1:N_req)

        # Run offline algorithm tests
        capacity = 5
        min_cost, end_t, opt_route = OfflineAlgorithm.run(g, requests, capacity)
        @test opt_route == [sort(request_vertices)]

        capacity = 1
        min_cost, end_t, opt_route = OfflineAlgorithm.run(g, requests, capacity)
        @test opt_route == map(el -> [el], sort(request_vertices))
    end
end

@testset "Online algorithm Wait And Ignore tests" begin
    @testset "Simple small graph test" begin
        # Create simple weighted graph
        srcs = [1, 1, 2]
        dsts = [2, 3, 3]
        wgts = [10., 4., 7.]
        g = SimpleWeightedGraph(srcs, dsts, wgts);

        # Define requests
        N_req = 2
        release_times = [0, 5]
        request_vertices = [2, 3]
        requests = map(i -> DataModel.Request(release_times[i], request_vertices[i]), 1:N_req)

        # Define capacity
        capacity = 2

        # Run Wait And Ignore algorithm tests
        min_cost, end_t, opt_route = OfflineAlgorithm.run(g, requests, capacity)
        cost, end_t, route = WaitAndIgnore.run(g, requests, capacity)

        @test cost ≈ 32.0
        @test end_t ≈ 33.0
        @test route == [[3], [2]]
        @test cost / min_cost ≈ 1.28
    end
end

@testset "Online algorithm Wait And Return tests" begin
    @testset "Simple small graph test" begin
        # Create simple weighted graph
        srcs = [1,   1,   1,   1,   1,   2,   3]
        dsts = [2,   3,   4,   5,   6,   3,   4]
        wgts = [10., 10., 40., 10., 10., 10., 20.]
        g = SimpleWeightedGraph(srcs, dsts, wgts);

        # Define requests
        N_req = 5
        release_times = [30, 30, 30, 45, 55]
        request_vertices = [2, 3, 4, 5, 6]
        requests = map(i -> DataModel.Request(release_times[i], request_vertices[i]), 1:N_req)

        # Define capacity
        capacity = 10
        
        # Run Wait And Return algorithm tests
        cost, end_t, route = WaitAndReturn.run(g, requests, capacity)

        @test cost ≈ 410.0
        @test end_t ≈ 170.0
        @test route == [[2, 3], [5, 6, 4]]
    end
end

@testset "MILP model for offline algorithm" begin
    @testset "Small graph with only one route" begin 
        # Create simple weighted graph
        srcs = [1,   1,   1,   2,   4 ]
        dsts = [2,   3,   4,   3,   3 ]
        wgts = [10., 20., 5.,  8.,  5.]
        g = SimpleWeightedGraph(srcs, dsts, wgts);

        # Define requests
        N_req = 2
        release_times = zeros(N_req)
        request_vertices = [2, 3]
        requests = map(i -> DataModel.Request(release_times[i], request_vertices[i]), 1:N_req)

        # Define capacity
        capacity = N_req

        # Run MILP offline algorithm tests
        min_cost, end_t = MilpOfflineAlgorithm.run(g, requests, capacity)
        @test min_cost ≈ 28.0
        @test end_t ≈ 28.0
    end

    @testset "Different capacities and different routes" begin
        # Create simple weighted graph
        srcs = [1,   1,   1,   1,   1,   2,   3]
        dsts = [2,   3,   4,   5,   6,   3,   4]
        wgts = [10., 10., 40., 10., 10., 10., 20.]
        g = SimpleWeightedGraph(srcs, dsts, wgts);

        # Define requests
        N_req = 3
        release_times = [0, 0, 0]
        request_vertices = [2, 3, 5]
        requests = map(i -> DataModel.Request(release_times[i], request_vertices[i]), 1:N_req)

        # Define capacity and run MILP offline algorithm tests
        capacity = 3
        min_cost, end_t = MilpOfflineAlgorithm.run(g, requests, capacity)
        @test min_cost ≈ 70.0
        @test end_t ≈ 50.0

        capacity = 2
        min_cost, end_t = MilpOfflineAlgorithm.run(g, requests, capacity)
        @test min_cost ≈ 70.0
        @test end_t ≈ 50.0

        capacity = 1
        min_cost, end_t = MilpOfflineAlgorithm.run(g, requests, capacity)
        @test min_cost ≈ 90.0
        @test end_t ≈ 60.0
    end

    @testset "Reverse path" begin
        # Create simple weighted graph
        srcs = [1,   1,     1,  2  ]
        dsts = [2,   3,     4,  3  ]
        wgts = [10., 1000., 100., 20.]
        g = SimpleWeightedGraph(srcs, dsts, wgts);

        # Define requests
        N_req = 2
        release_times = zeros(N_req)
        request_vertices = [3, 4]
        requests = map(i -> DataModel.Request(release_times[i], request_vertices[i]), 1:N_req)

        # Define capacity
        capacity = 1

        # Run MILP offline algorithm tests
        min_cost, end_t = MilpOfflineAlgorithm.run(g, requests, capacity)

        @test min_cost ≈ 190.0
        @test end_t ≈ 260.0
    end

    @testset "Initial time" begin
        # Create simple weighted graph
        srcs = [1,   1,     1,  2  ]
        dsts = [2,   3,     4,  3  ]
        wgts = [10., 1000., 100., 20.]
        g = SimpleWeightedGraph(srcs, dsts, wgts);

        # Define requests
        N_req = 2
        release_times = zeros(N_req)
        request_vertices = [3, 4]
        requests = map(i -> DataModel.Request(release_times[i], request_vertices[i]), 1:N_req)

        # Define capacity
        capacity = 1

        # Define initial time
        initial_t = 50.0

        # Run MILP offline algorithm tests
        min_cost, end_t = MilpOfflineAlgorithm.run(g, requests, capacity, initial_t)

        @test min_cost ≈ 290.0
        @test end_t ≈ 310.0
    end

    @testset "Random euclidean graph instance" begin
        # Create an euclidean graph 
        g, dists = euclidean_graph(6, 2, p=2)
        srcs = Vector{Int}()
        dsts = Vector{Int}()
        weights = Vector{Float64}()
        for e in dists
            push!(srcs, e.first.src)
            push!(dsts, e.first.dst)
            push!(weights, e.second)
        end
        g = SimpleWeightedGraph(srcs, dsts, weights);
        
        # Define requests
        N_req = 4
        max_weight = maximum(weights)
        release_times = repeat([max_weight], N_req)
        request_vertices = sample(2:nv(g), N_req ,replace=false)
        requests = map(i -> DataModel.Request(release_times[i], request_vertices[i]), 1:N_req)

        # Define initial time
        initial_t = max_weight

        # Test for different capacities
        for capacity in 1:N_req
            # Run MILP offline algorithm tests
            min_cost_milp, end_t_milp = MilpOfflineAlgorithm.run(g, requests, capacity, initial_t)
            min_cost_brute_force, end_t_brute, opt_route = OfflineAlgorithm.run(g, requests, capacity, initial_t)

            @test min_cost_milp ≈ min_cost_brute_force
        end
    end
end