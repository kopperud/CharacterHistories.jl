export Node
export Tip
export Branch
export Root
export tiplabels

abstract type AbstractNode end
abstract type AbstractBranch end
abstract type InternalNode <: AbstractNode end

mutable struct Tip <: AbstractNode
    index::Int64
    inbounds::AbstractBranch
    species_name::String
    Tip() = new()
end

mutable struct Branch <: AbstractBranch
    index::Int64
    inbounds::AbstractNode
    outbounds::AbstractNode

    states::Vector{String}
    times::Vector{Float64}
    Branch() = new()
end

mutable struct Node <: InternalNode
    index::Int64
    inbounds::AbstractBranch
    left::AbstractBranch
    right::AbstractBranch
    Node() = new()
end

mutable struct Root <: InternalNode
    index::Int64
    left::AbstractBranch
    right::AbstractBranch
    Root() = new()
end



function tl_postorder!(labels::Vector{String}, node::Tip)
    index::Int64
    label = node.species_name
    push!(labels, label)
end


function tiplabels(node::Root)
    labels = String[]

    left_branch = node.left
    right_branch = node.right

    left_node = left_branch.outbounds
    right_node = right_branch.outbounds

    tl_postorder!(labels, left_node)
    tl_postorder!(labels, right_node)
end

function tl_postorder!(labels::Vector{String}, node::InternalNode)
    left_branch = node.left
    right_branch = node.right

    left_node = left_branch.outbounds
    right_node = right_branch.outbounds

    tl_postorder!(labels, left_node)
    tl_postorder!(labels, right_node)
end

export tipstates

function tipstates(tree::Root)
    data = Dict{String,String}()

    data = ts_postorder!(tree, data)
end

function ts_postorder!(node::T, data::Dict{String,String}) where {T <: InternalNode}
    left_branch = node.left
    right_branch = node.right

    left_node = left_branch.outbounds
    right_node = right_branch.outbounds

    data = ts_postorder!(left_node, data)
    data = ts_postorder!(right_node, data)
end

function ts_postorder!(node::Tip, data::Dict{String,String})
    label = node.species_name
    parent_branch = node.inbounds
    most_recent_state = parent_branch.states[1]

    data[label] = most_recent_state
    return(data)
end