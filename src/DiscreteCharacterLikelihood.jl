export likelihood, postorder

function likelihood(tree, model, data)
    n_branches = number_of_edges(tree)
    D = zeros(n_branches, 2, model.k)

    x, log_nf = postorder!(tree, model, data, D)

    root_freqs = [1/model.k for _ in 1:model.k]

    

    logl = log(sum(root_freqs .* x)) + log_nf
    return(logl, D)
end

function postorder!(node::T, model::Mk, data::Dict{String,Int64}, D::Array{Float64, 3}) where {T<:InternalNode}
    left_branch = node.left
    left_node = left_branch.outbounds
    left_bl = sum(left_branch.times)
    P_left = transition_probability(model, left_bl)

    right_branch = node.right
    right_node = right_branch.outbounds
    right_bl = sum(right_branch.times)
    P_right = transition_probability(model, right_bl)

    x_left, log_nf_left = postorder!(left_node, model, data, D)
    x_right, log_nf_right = postorder!(right_node, model, data, D)


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

function postorder!(node::Tip, model::Mk, data::Dict, D::Array{Float64, 3})
    state = zeros(model.k)
    state[data[node.species_name]] = 1
    log_nf = 0.0

    return(state, log_nf)
end

export ancestral_state_probabilities

function ancestral_state_probabilities(tree, model, data, D)
    n_branches = number_of_edges(tree)
    F = zeros(n_branches, 2, model.k)

    left_branch_idx = tree.left.index
    right_branch_idx = tree.right.index

    F0 = ones(model.k)
    D0 = D[left_branch_idx,2,:] .* D[right_branch_idx,2,:]
    S0 = F0 .* D0
    S0 = S0 ./ sum(S0)

    preorder!(tree, model, data, D, F, S0)

    
    S = D[:,1,:] .* F[:,1,:]
    return(S)
#    x, log_nf = postorder!(tree, model, data, D)

#    root_freqs = [1/model.k for _ in 1:model.k]

#    logl = log(sum(root_freqs .* x)) + log_nf
#    return(logl, D)
end

function preorder!(node::T, model::Mk, data::Dict{String,Int64}, D::Array{Float64, 3}, F::Array{Float64, 3}, S_parent) where {T <: InternalNode}
    left_branch = node.left
    left_node = left_branch.outbounds
    left_bl = sum(left_branch.times)
    P_left = transition_probability(model, left_bl)
    F_left = S_parent ./ D[left_branch.index,2,:]
    F[left_branch.index,2,:] = F_left
    F[left_branch.index,1,:] = P_left * F_left
    S_left = F[left_branch.index,1,:] .* D[left_branch.index,1,:]
    S_left = S_left ./ sum(S_left)

    right_branch = node.right
    right_node = right_branch.outbounds
    right_bl = sum(right_branch.times)
    P_right = transition_probability(model, right_bl)
    F_right = S_parent ./ D[right_branch.index,2,:]
    F[right_branch.index,2,:] = F_right
    F[right_branch.index,1,:] = P_right * F_right
    S_right = F[right_branch.index,1,:] .* D[right_branch.index,1,:]
    S_right = S_right ./ sum(S_right)

    preorder!(left_node, model, data, D, F, S_left)
    preorder!(right_node, model, data, D, F, S_right)
end

function preorder!(node::Tip, model::Mk, data::Dict, D::Array{Float64, 3}, F::Array{Float64, 3}, S::Vector{Float64})
end