export redraw_nodes_recursive!

function redraw_nodes_recursive!(
        node::T, 
        model::CharacterHistory, 
        S::Array{Float64, 2}
    ) where {T <: InternalNode}
    
    ## redraw this node
    d = Distributions.Categorical(S[node.index,:])
    r = rand(d)

    new_state = model.state_space[r]
    node.state = new_state

    left_node = node.left.outbounds
    redraw_nodes_recursive!(left_node, model, S)

    right_node = node.right.outbounds
    redraw_nodes_recursive!(right_node, model, S)
end

function redraw_nodes_recursive!(node::Tip, model::CharacterHistory, S::Array{Float64, 2})
    ## dont do anything
end


function redraw_branches_recursive!(node::Root, model::CharacterHistory)
    oldest_state = node.state

    ## left branch
    left_branch = node.left
    left_node = left_branch.outbounds
    left_state = left_node.state
    redraw_branch!(left_branch, model, oldest_state, left_state)

    ## right branch
    right_branch = node.right
    right_node = right_branch.outbounds
    right_state = right_node.state
    redraw_branch!(right_branch, model, oldest_state, right_state)

    redraw_branches_recursive!(left_node, model)
    redraw_branches_recursive!(right_node, model)
end

function redraw_branches_recursive!(node::Node, model::CharacterHistory)
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


    redraw_branches_recursive!(left_node, model)
    redraw_branches_recursive!(right_node, model)
end

function redraw_branches_recursive!(node::Tip, model::CharacterHistory) 
    ## dont do anything
end


