export brownian_stateless

function brownian_stateless(tree::Root, model::Brownian, data::Dict)
    n_nodes = number_of_nodes(tree)

    μ = zeros(n_nodes)
    V = zeros(n_nodes)


    log_nf_factor = brownian_stateless_po(tree, model, μ, V, data)

    ## root treatment
    ## assume that x_0 = μ_root
    root_index = tree.index
    #log_likelihood = log_nf_factor - 0.5 * log(pi * V[root_index])
    d = Distributions.Normal(μ[root_index], sqrt(V[root_index]))
    log_likelihood = log_nf_factor + Distributions.logpdf(d, model.mean)
    return(log_likelihood)#, μ, V)
end

function brownian_stateless_po(
        node::T,
        model::Brownian, 
        μ::Vector{Float64},
        V::Vector{Float64}, 
        data::Dict
    ) where {T <: InternalNode}

    left_branch = node.left
    left_node = left_branch.outbounds
    left_bl = sum(left_branch.times)

    right_branch = node.right
    right_node = right_branch.outbounds
    right_bl = sum(right_branch.times)

    log_nf_left = brownian_stateless_po(left_node, model, μ, V, data)
    log_nf_right = brownian_stateless_po(right_node, model, μ, V, data)

    μ_left = μ[left_node.index]
    V_left = V[left_node.index] + left_bl * model.sigma2

    μ_right = μ[right_node.index]
    V_right = V[right_node.index] + right_bl * model.sigma2

    ## merging rule
    μ_node = (μ_left*V_right + μ_right*V_left) / (V_left+V_right)
    V_node = (V_left*V_right)/(V_left+V_right)

    #println(μ_left)
    μ[node.index] = μ_node
    V[node.index] = V_node

    ## normalizing factor
    contrast = μ_left - μ_right
    log_nf = -(1.0)*contrast^2 / (2*(V_left+V_right))
    log_nf -= 0.5 * log(2 * pi * (V_left+V_right))
    log_nf += log_nf_left
    log_nf += log_nf_right

    return(log_nf)
end

function brownian_stateless_po(node::Tip, model::Brownian, μ::Vector{Float64}, V::Vector{Float64}, data::Dict)
    tip_label = node.species_name

    μ[node.index] = data[tip_label]
    V[node.index] = 0.0
    #V = 0.0 ## assume no measurement error
    log_nf = 0.0

    return(log_nf)
end
