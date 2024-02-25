export loglikelihood

function loglikelihood(
    tree::Root, 
    model::BrownianSD, 
    data::Dict{String,T}
    ) where {T <: Real}
    n_nodes = number_of_nodes(tree)

    μ = zeros(eltype(model.mean), n_nodes)
    V = zeros(eltype(model.mean), n_nodes)

    log_nf_factor = loglikelihood_po(tree, model, μ, V, data)

    ## root treatment
    ## assume that x_0 = μ_root
    root_index = tree.index
    #root_state = tree.left.states[end]

    d = Distributions.Normal(μ[root_index], sqrt(V[root_index]))
    log_likelihood = log_nf_factor + Distributions.logpdf(d, model.mean)
    return(log_likelihood)
end

function loglikelihood_po(
        node::N,
        model::BrownianSD, 
        μ::Vector{T},
        V::Vector{T},
        data::Dict
    ) where {N <: InternalNode, T <: Real}

    left_branch = node.left
    left_node = left_branch.outbounds

    right_branch = node.right
    right_node = right_branch.outbounds

    log_nf_left = loglikelihood_po(left_node, model, μ, V, data)
    log_nf_right = loglikelihood_po(right_node, model, μ, V, data)

    μ_left = μ[left_node.index]
    V_left = V[left_node.index]
    for (state, time) in zip(left_branch.states, left_branch.times)
        state_idx = findfirst(isequal(state), model.state_space)
        V_left += time * model.sigma2[state_idx]
    end    

    μ_right = μ[right_node.index]
    V_right = V[right_node.index]
    for (state, time) in zip(right_branch.states, right_branch.times)
        state_idx = findfirst(isequal(state), model.state_space)
        V_right += time * model.sigma2[state_idx]
    end

    ## merging rule
    μ_node = (μ_left*V_right + μ_right*V_left) / (V_left+V_right)
    V_node = (V_left*V_right)/(V_left+V_right)

    μ[node.index] = μ_node
    V[node.index] = V_node

    ## normalizing factor
    contrast = μ_left - μ_right
    log_nf = -(1.0)*contrast^2 / (2*(V_left+V_right))
    PI = 3.14159265358979
    log_nf -= 0.5 * log(2.0 * PI * (V_left+V_right))
    log_nf += log_nf_left
    log_nf += log_nf_right

    return(log_nf)
end

function loglikelihood_po(
    node::Tip, 
    model::BrownianSD, 
    μ::Vector{T}, 
    V::Vector{T},
    data::Dict
    )where {T <: Real}

    tip_label = node.species_name

    μ[node.index] = data[tip_label]
    V[node.index] = zero(T)
    #V = 0.0 ## assume no measurement error
    log_nf = zero(T)

    return(log_nf)
end
