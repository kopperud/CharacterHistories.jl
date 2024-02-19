export likelihood, postorder

function likelihood(tree, model, data)
    x, log_nf = postorder(tree, model, data)

    root_freqs = [1/model.k for _ in 1:model.k]

    logl = log(sum(root_freqs .* x)) + log_nf
    return(logl)
end

function postorder(node::T, model::Mk, data::Dict{String,Int64}) where {T<:InternalNode}
    left_branch = node.left
    left_node = left_branch.outbounds
    left_bl = sum(left_branch.times)
    P_left = transition_probability(model, left_bl)

    right_branch = node.right
    right_node = right_branch.outbounds
    right_bl = sum(right_branch.times)
    P_right = transition_probability(model, right_bl)

    x_left, log_nf_left = postorder(left_node, model, data)
    x_right, log_nf_right = postorder(right_node, model, data)


    state = P_left * x_left .* P_right * x_right
    nf_node = sum(state)
    state = state ./ nf_node
    log_nf = log_nf_left + log_nf_right + log(nf_node)

    return(state, log_nf)
end

function postorder(node::Tip, model::Mk, data::Dict)
    state = zeros(model.k)
    state[data[node.species_name]] = 1
    log_nf = 0.0

    return(state, log_nf)
end
