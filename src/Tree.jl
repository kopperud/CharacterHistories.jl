export Node
export Tip
export Branch
export Root

abstract type AbstractNode end
abstract type AbstractBranch end

mutable struct Tip
    node_index::Int64
    species_names::String
end

mutable struct Branch <: AbstractBranch
    inbounds::Union{Nothing,AbstractNode}
    outbounds::Union{Nothing,AbstractNode,Tip}

    states::Vector{String}
    times::Vector{Float64}
end

mutable struct Node <: AbstractNode
    inbounds::Union{Nothing,AbstractBranch}
    left::Union{Nothing,AbstractBranch,Tip}
    right::Union{Nothing,AbstractBranch,Tip}
end

mutable struct Root <: AbstractNode
    left::Union{Nothing,AbstractBranch,Tip}
    right::Union{Nothing,AbstractBranch,Tip}
end





