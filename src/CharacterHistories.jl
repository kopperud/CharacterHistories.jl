module CharacterHistories

import ForwardDiff
import Distributions
import LinearAlgebra
import ProgressMeter
import Base
import OrderedCollections
import Random

greet() = print("Hello World!")

include("tree/types.jl")
include("tree/Tree.jl")
include("tree/reindex.jl")
include("tree/number_of_nodes.jl")
include("tree/copy.jl")
include("tree/number_state_changes.jl")
include("tree/treelength.jl")
include("tree/show.jl")

abstract type AbstractMove end

include("random_variables/random_variables.jl")
include("random_variables/dag.jl")
include("random_variables/normal.jl")
include("random_variables/constant.jl")
include("random_variables/exponential.jl")


include("discrete_character/Mk_model.jl")
include("continuous_character/models.jl")



include("io/readnewick.jl")
include("io/writenewick.jl")

include("discrete_character/preorder.jl")
include("discrete_character/postorder.jl")
include("discrete_character/likelihood.jl")
include("discrete_character/ancestral_state_probabilities.jl")

include("proposals/stochastic_character_map.jl")
include("proposals/sample_branch.jl")
include("proposals/redraw_recursive.jl")
include("proposals/redraw_node.jl")
include("proposals/redraw_branch.jl")
include("proposals/reassign_copies.jl")

include("moves/moves.jl")
include("moves/scale.jl")
include("moves/slide.jl")


include("mcmc/mcmc.jl")
include("mcmc/ess.jl")

include("continuous_character/brownian_stateless.jl")
include("continuous_character/brownian_state_dependent.jl")
include("continuous_character/ou_stateless.jl")
include("continuous_character/ou_state_dependent.jl")
include("continuous_character/vcv.jl")

include("find_node.jl")
include("parent_nodes.jl")

include("gradient_descent.jl")

include("plot.jl")
include("utils.jl")

# Path into package
export path
path(x...; dir::String = "data") = joinpath(@__DIR__, "..", dir, x...)

end # module CharacterHistories
