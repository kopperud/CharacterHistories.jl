export Node
export Tip
export Branch
export Root
export InternalNode

abstract type AbstractNode end
abstract type AbstractBranch end
abstract type InternalNode <: AbstractNode end

mutable struct Tip <: AbstractNode
    index::Int64
    inbounds::AbstractBranch
    species_name::String
    state::String
    Tip() = new()
end

mutable struct Branch <: AbstractBranch
    index::Int64
    inbounds::AbstractNode
    outbounds::AbstractNode

    ## the [1]   first item is the most recent state
    ## the [end] final item is the oldest state
    states::Vector{String}
    times::Vector{Float64}
    Branch() = new()
end

mutable struct Node <: InternalNode
    index::Int64
    inbounds::AbstractBranch
    left::AbstractBranch
    right::AbstractBranch
    state::String
    Node() = new()
end

mutable struct Root <: InternalNode
    index::Int64
    left::AbstractBranch
    right::AbstractBranch
    state::String
    Root() = new()
end