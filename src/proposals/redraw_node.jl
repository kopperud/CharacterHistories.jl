
export redraw_node!

function redraw_node!(node::Node, model::Mk)
    this_node_state = node.state

    ## parent branch
    parent_branch = node.inbounds
    parent_node = parent_branch.inbounds
    parent_state = parent_node.state
    redraw_branch!(parent_branch, model, parent_state, this_node_state)

    ## left branch
    left_branch = node.left
    left_node = left_branch.outbounds
    left_state = left_node.state
    redraw_branch!(left_branch, model, this_node_state, left_state)

    ## right branch
    right_branch = node.right
    right_node = right_branch.outbounds
    right_state = right_node.state
    redraw_branch!(right_branch, model, this_node_state, right_state)
end

function redraw_node!(node::Root, model::Mk)
    this_node_state = node.state

    ## left branch
    left_branch = node.left
    left_node = left_branch.outbounds
    left_state = left_node.state
    redraw_branch!(left_branch, model, this_node_state, left_state)

    ## right branch
    right_branch = node.right
    right_node = right_branch.outbounds
    right_state = right_node.state
    redraw_branch!(right_branch, model, this_node_state, right_state)
end