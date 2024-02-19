export Node
export Tip
export Branch
export Root

abstract type AbstractNode end
abstract type AbstractBranch end
abstract type InternalNode <: AbstractNode end

mutable struct Tip <: AbstractNode
    node_index::Int64
    species_name::String
end

mutable struct Branch <: AbstractBranch
    inbounds::Union{Nothing,AbstractNode} ## this union with Nothing cause type unstability/branches, not good
    outbounds::Union{Nothing,AbstractNode} ## how else to initialize the branch if I dont have the parent?

    states::Vector{String}
    times::Vector{Float64}
end

mutable struct Node <: InternalNode
    inbounds::Union{Nothing,AbstractBranch}
    left::Union{Nothing,AbstractBranch}
    right::Union{Nothing,AbstractBranch}
end

mutable struct Root <: InternalNode
    left::Union{Nothing,AbstractBranch}
    right::Union{Nothing,AbstractBranch}
end





