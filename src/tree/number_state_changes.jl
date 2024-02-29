
export number_of_state_changes

function number_of_state_changes(tree::Root)
    n = Int64[0]
     
    number_of_state_changes!(tree, n)

    return(n[1])
end

function number_of_state_changes!(
        node::T,
        n::Vector{Int64}
    ) where {T <: InternalNode}

    left_branch = node.left
    left_node = left_branch.outbounds
    n_left_changes = length(left_branch.times) - 1
    
    n[1] += n_left_changes

    right_branch = node.right
    right_node = right_branch.outbounds
    n_right_changes = length(right_branch.times) - 1
    
    n[1] += n_right_changes

    number_of_state_changes!(left_node, n)
    number_of_state_changes!(right_node, n)    
end

function number_of_state_changes!(
    node::Tip,
    n::Vector{Int64}
) 
    ## do nothing
end