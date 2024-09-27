
export redraw_node!

function redraw_node!(
        node::Node, 
        model::CharacterHistory,
        S::Array{Float64, 2}
    )

    ## should the node state be drawn uniformly?
    d = Distributions.Categorical(S[node.index,:])
    r = rand(d)

    new_state = model.state_space[r]
    node.state = new_state

    ## parent branch
    parent_branch = node.inbounds
    parent_node = parent_branch.inbounds
    parent_state = parent_node.state
    redraw_branch!(parent_branch, model, parent_state, new_state)

    ## left branch
    left_branch = node.left
    left_node = left_branch.outbounds
    left_state = left_node.state
    redraw_branch!(left_branch, model, new_state, left_state)

    ## right branch
    right_branch = node.right
    right_node = right_branch.outbounds
    right_state = right_node.state
    redraw_branch!(right_branch, model, new_state, right_state)
end

function redraw_node!(
        node::Root, 
        model::CharacterHistory, 
        S::Array{Float64, 2}
    )
    #this_node_state = node.state
    d = Distributions.Categorical(S[node.index,:])
    r = rand(d)

    new_state = model.state_space[r]
    node.state = new_state

    ## left branch
    left_branch = node.left
    left_node = left_branch.outbounds
    left_state = left_node.state
    redraw_branch!(left_branch, model, new_state, left_state)

    ## right branch
    right_branch = node.right
    right_node = right_branch.outbounds
    right_state = right_node.state
    redraw_branch!(right_branch, model, new_state, right_state)
end