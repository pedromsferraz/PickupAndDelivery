module FinalProjectPedroFerraz

export DataModel, OfflineAlgorithm

include("DataModel.jl")
using .DataModel

include("OfflineAlgorithm.jl")
using .OfflineAlgorithm

end # module
