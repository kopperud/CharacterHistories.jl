module CharacterHistories

greet() = print("Hello World!")

include("Tree.jl")
include("Newick.jl")

# Path into package
export path
path(x...; dir::String = "data") = joinpath(@__DIR__, "..", dir, x...)

end # module CharacterHistories
