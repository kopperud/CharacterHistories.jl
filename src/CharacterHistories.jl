module CharacterHistories

import ForwardDiff
import Distributions
import LinearAlgebra
import ProgressMeter
import Base
import OrderedCollections
import Random

greet() = print("Hello World!")

## the tree data structure
include("tree/types.jl")
include("tree/Tree.jl")
include("tree/reindex.jl")
include("tree/number_of_nodes.jl")
include("tree/copy.jl")
include("tree/number_state_changes.jl")
include("tree/treelength.jl")
include("tree/show.jl")
include("tree/find_node.jl")
include("tree/parent_nodes.jl")

## dag
include("dag/abstract_types.jl")
include("dag/dag.jl")
include("dag/constant.jl")



## distributions
include("distributions/exponential.jl")
include("distributions/normal.jl")
include("distributions/lognormal.jl")

include("distributions/discrete_character/phyloctmc_dataaugmented.jl")

include("distributions/continuous_character/brownian_state_dependent.jl")
include("distributions/continuous_character/brownian_stateless.jl")
include("distributions/continuous_character/ou_stateless.jl")
include("distributions/continuous_character/ou_state_dependent.jl")
include("distributions/continuous_character/vcv.jl")

## moves
include("moves/moves.jl")
include("moves/slide.jl")
include("moves/scale.jl")

## tree moves
include("moves/characterhistory/sample_branch.jl")
include("moves/characterhistory/redraw_recursive.jl")
include("moves/characterhistory/redraw_node.jl")
include("moves/characterhistory/redraw_branch.jl")
include("moves/characterhistory/reassign_copies.jl")

## draw scm
include("proposals/stochastic_character_map.jl")

## write and read trees
include("io/readnewick.jl")
include("io/writenewick.jl")

## mcmc 
include("mcmc/mcmc.jl")
include("mcmc/ess.jl")

## gd
include("gradient_descent.jl")

## extras
include("plot.jl")
include("utils.jl")

# Path into package
export path
path(x...; dir::String = "data") = joinpath(@__DIR__, "..", dir, x...)

end # module CharacterHistories
