export tiplabels




function tl_postorder!(labels::Vector{String}, node::Tip)
    #index::Int64
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

function get_index(tree::Root)
    return(tree.index)
end

function get_index(node::Node)
    return(node.index)
end

function get_index(tip::Tip)
    return(tip.index)
end

function get_index(branch::Branch)
    return(branch.index)
end

export get_state_space

function get_state_space(tree::Root)
    s = Set{String}()

    state_space!(tree, s)
    v = collect(s)
    return(v)
end

function state_space!(node::N, s::Set{String}) where {N <: InternalNode}
    left_branch = node.left
    right_branch = node.right

    for state in left_branch.states
        push!(s, state)
    end
    
    for state in right_branch.states
        push!(s, state)
    end

    state_space!(left_branch.outbounds, s)
    state_space!(right_branch.outbounds, s)
end

function state_space!(node::Tip, s::Set{String})
    ## do nothing    
end