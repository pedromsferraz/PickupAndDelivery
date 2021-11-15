using Graphs, SimpleWeightedGraphs
using GraphPlot, Cairo, Compose

# Create a random weighted graph with N vertices and M edges
N = 5
M = 7
graph = SimpleDiGraph(N, M);
g = SimpleWeightedDiGraph(N);
weights = Vector{Float64}()
for e in edges(graph)
    w = rand((1:20))
    add_edge!(g, e.src, e.dst, w)
    push!(weights, w);
end

using GraphPlot, Cairo, Compose
# Plot the graph
mkpath("images")
plot = gplot(g, nodelabel=1:nv(g), 
                edgelabel=weights);
draw(PDF("images/graph.pdf", 16cm, 16cm), plot);

# Create an euclidean graph
g, dists = euclidean_graph(10, 2, p=2)
srcs = Vector{Int}()
dsts = Vector{Int}()
weights = Vector{Float64}()
for e in dists
    push!(srcs, e.first.src)
    push!(dsts, e.first.dst)
    push!(weights, e.second)
end
g = SimpleWeightedGraph(srcs, dsts, weights);

plot = gplot(g, nodelabel=1:nv(g), 
                edgelabel=weights);
draw(PDF("images/graph.pdf", 32cm, 32cm), plot);

# Calculate and plot MST using Kruskal's algorithm
mst_edges = kruskal_mst(g)
srcs = map(e -> e.src, mst_edges)
dsts = map(e -> e.dst, mst_edges)
weights = map(e -> e.weight, mst_edges)
mst = SimpleWeightedGraph(srcs, dsts, weights)

plot = gplot(mst, nodelabel=1:nv(mst), 
                edgelabel=weights);
draw(PDF("images/MST.pdf", 32cm, 32cm), plot);
