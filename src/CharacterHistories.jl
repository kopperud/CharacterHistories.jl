module CharacterHistories

import Distributions

greet() = print("Hello World!")

include("Tree.jl")
include("io/readnewick.jl")
include("io/writenewick.jl")

include("discrete_character/preorder.jl")
include("discrete_character/postorder.jl")
include("discrete_character/likelihood.jl")
include("discrete_character/ancestral_state_probabilities.jl")
include("discrete_character/Mk_model.jl")

include("reindex.jl")
include("number_of_nodes.jl")
include("plot.jl")
include("StochasticCharacterMap.jl")

# Path into package
export path
path(x...; dir::String = "data") = joinpath(@__DIR__, "..", dir, x...)

end # module CharacterHistories
