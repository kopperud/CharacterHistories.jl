export tree_length

function tree_length(node::N) where {N <: InternalNode}
    left_node = node.left.outbounds
    right_node = node.right.outbounds

    tl_left = tree_length(left_node) + sum(node.left.times)
    tl_right = tree_length(right_node) + sum(node.right.times)

    return(tl_left+tl_right)
end

function tree_length(node::Tip)
    return(0.0)
end
