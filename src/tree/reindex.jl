
export reindex!

function reindex!(tree::Root)
    tip_index = Int64[1]
    branch_index = Int64[1]
    
    ## first index the tip nodes, and return the ntip+1 index
    root_index = tipindex_po!(tree, tip_index)

    ## secondly, index the branches and the internal nodes, 
    ## starting from the root and going (left, right)
    branchindex_po!(tree, root_index, branch_index)
end

function tipindex_po!(
    node::T, 
    tip_index::Vector{Int64}
    ) where {T<:InternalNode}

    left_branch = node.left
    left_node = left_branch.outbounds

    tipindex_po!(left_node, tip_index)

    right_branch = node.right
    right_node = right_branch.outbounds

    tipindex_po!(right_node, tip_index)

    return(tip_index)
end

function tipindex_po!(
    node::Tip, 
    tip_index::Vector{Int64}
    )

    node.index = tip_index[1]
    tip_index[1] += 1
end


function branchindex_po!(
    node::T, 
    node_index::Vector{Int64}, 
    branch_index::Vector{Int64}
    ) where {T<:InternalNode}

    node.index = node_index[1]
    node_index[1] += 1

    left_branch = node.left
    left_branch.index = branch_index[1]
    branch_index[1] += 1
    left_node = left_branch.outbounds

    branchindex_po!(left_node, node_index, branch_index)

    right_branch = node.right
    right_branch.index = branch_index[1]
    branch_index[1] += 1
    right_node = right_branch.outbounds

    branchindex_po!(right_node, node_index, branch_index)


end

function branchindex_po!(::Tip, ::Vector{Int64}, ::Vector{Int64})
end
