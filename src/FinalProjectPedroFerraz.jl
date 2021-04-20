module FinalProjectPedroFerraz

export DataModel, OfflineAlgorithm, Graph

include("DataModel.jl")
using .DataModel

include("Graph.jl")
using .Graph

include("OfflineAlgorithm.jl")
using .OfflineAlgorithm

end # module
