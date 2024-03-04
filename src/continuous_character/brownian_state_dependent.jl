export loglikelihood

const PI = 3.14159265358979

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
    root_index = get_index(tree)
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

    left_index = get_index(left_node)
    μ_left = μ[left_index]
    V_left = V[left_index]
    for (time, state) in zip(left_branch.times, left_branch.states)
        sigma2 = model.sigma2[state]
        V_left += time * sigma2 ## not type stable
    end

    μ_right = μ[get_index(right_node)]
    V_right = V[get_index(right_node)]
    for (time, state) in zip(right_branch.times, right_branch.states)
        sigma2 = model.sigma2[state]
        V_right += time * sigma2 ## not type stable
    end

    ## merging rule
    μ_node = (μ_left*V_right + μ_right*V_left) / (V_left+V_right)
    V_node = (V_left*V_right)/(V_left+V_right)

    node_index = get_index(node)
    μ[node_index] = μ_node
    V[node_index] = V_node

    ## normalizing factor
    contrast = μ_left - μ_right
    log_nf = -(1.0)*contrast^2 / (2*(V_left+V_right))
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

    node_index = get_index(node)
    μ[node_index] = data[tip_label]
    V[node_index] = zero(T)
    #V = 0.0 ## assume no measurement error
    log_nf = zero(T)

    return(log_nf)
end
