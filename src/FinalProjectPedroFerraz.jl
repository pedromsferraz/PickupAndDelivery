module FinalProjectPedroFerraz

export DataModel, GraphPreprocessing, OfflineAlgorithm, WaitAndIgnore, WaitAndReturn

include("DataModel.jl")
using .DataModel

include("GraphPreprocessing.jl")
using .GraphPreprocessing

include("OfflineAlgorithm.jl")
using .OfflineAlgorithm

include("WaitAndIgnore.jl")
using .WaitAndIgnore

include("WaitAndReturn.jl")
using .WaitAndReturn

end # module
