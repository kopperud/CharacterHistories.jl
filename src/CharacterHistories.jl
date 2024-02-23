module CharacterHistories

import Distributions

greet() = print("Hello World!")

include("Tree.jl")
include("Newick.jl")
include("SubstitutionModel.jl")
include("DiscreteCharacterLikelihood.jl")
include("reindex.jl")
include("number_of_nodes.jl")
include("plot.jl")
include("StochasticCharacterMap.jl")

# Path into package
export path
path(x...; dir::String = "data") = joinpath(@__DIR__, "..", dir, x...)

end # module CharacterHistories
