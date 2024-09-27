#######################################
##
##               postorder
##
########################################
function postorder!(
    node::N, 
    ch::CharacterHistory,
    D::Array{T, 3}
    ) where {N <: InternalNode, T <: Real}

    #α = getvalue(ch.α)

    left_branch = node.left
    left_node = left_branch.outbounds
    left_bl = sum(left_branch.times)
    P_left = transition_probability(ch, left_bl)

    right_branch = node.right
    right_node = right_branch.outbounds
    right_bl = sum(right_branch.times)
    P_right = transition_probability(ch, right_bl)

    x_left, log_nf_left = postorder!(left_node, ch, D)
    x_right, log_nf_right = postorder!(right_node, ch, D)

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
        ch::CharacterHistory,
        #model::CharacterHistory, 
        D::Array{T, 3}
    ) where {T <: Real}

    state = zeros(ch.k)
    idx = findfirst(isequal(node.state), ch.state_space)
    state[idx] = 1
    log_nf = 0.0

    return(state, log_nf)
end