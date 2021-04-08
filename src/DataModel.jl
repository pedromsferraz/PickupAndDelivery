module DataModel

import Base: isless, ==
export Request

mutable struct Request
    # Time at which request is released
    release_time::AbstractFloat

    # Destination vertex in graph
    destination::Int64
end

(==)(lhs::Request, rhs::Request) = (lhs.release_time == rhs.release_time) && (lhs.destination == rhs.destination)
isless(a::Request, b::Request) = isless(a.release_time, b.release_time)

end # module