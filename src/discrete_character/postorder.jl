#######################################
##
##               postorder
##
########################################
function postorder!(
    node::N, 
    model::Mk, 
    D::Array{T, 3}
    ) where {N <: InternalNode, T <: Real}
    left_branch = node.left
    left_node = left_branch.outbounds
    left_bl = sum(left_branch.times)
    P_left = transition_probability(model, left_bl)

    right_branch = node.right
    right_node = right_branch.outbounds
    right_bl = sum(right_branch.times)
    P_right = transition_probability(model, right_bl)

    x_left, log_nf_left = postorder!(left_node, model, D)
    x_right, log_nf_right = postorder!(right_node, model, D)

    D[left_branch.index,1,:] = x_left
    D[right_branch.index,1,:] = x_right
    
    state_left = P_left * x_left
    state_right = P_right * x_right

    D[left_branch.index,2,:] = state_left
    D[right_branch.index,2,:] = state_right
    
    state = state_left .* state_right
    nf_node = sum(state)
    state = state ./ nf_node
    log_nf = log_nf_left + log_nf_right + log(nf_node)

    return(state, log_nf)
end

function postorder!(
        node::Tip, 
        model::Mk, 
        D::Array{T, 3}
    ) where {T <: Real}

    state = zeros(model.k)
    idx = findfirst(isequal(node.state), model.state_space)
    state[idx] = 1
    log_nf = 0.0

    return(state, log_nf)
end