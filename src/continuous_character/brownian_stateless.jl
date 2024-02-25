export brownian
function brownian(
        tree::Root, 
        model::Brownian,
        data::Dict{String, T}
    ) where {T <: Real}
    n_nodes = number_of_nodes(tree)

    μ = zeros(T, n_nodes)
    V = zeros(T, n_nodes)

    log_nf_factor = brownian_po(tree, model, μ, V, data)

    ## root treatment
    ## assume that x_0 = μ_root
    root_index = tree.index
    #log_likelihood = log_nf_factor - 0.5 * log(pi * V[root_index])
    d = Distributions.Normal(μ[root_index], sqrt(V[root_index]))
    log_likelihood = log_nf_factor + Distributions.logpdf(d, model.mean)
    return(log_likelihood)#, μ, V)
end

function brownian_po(
        node::N,
        model::Brownian, 
        μ::Vector{T},
        V::Vector{T}, 
        data::Dict{String, T}
    ) where {N <: InternalNode, T <: Real}

    left_branch = node.left
    left_node = left_branch.outbounds
    left_bl = sum(left_branch.times)

    right_branch = node.right
    right_node = right_branch.outbounds
    right_bl = sum(right_branch.times)

    log_nf_left = brownian_po(left_node, model, μ, V, data)
    log_nf_right = brownian_po(right_node, model, μ, V, data)

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

function brownian_po(
    node::Tip, 
    model::Brownian, 
    μ::Vector{T},
    V::Vector{T}, 
    data::Dict{String, T}
    ) where {T <: Real}
    tip_label = node.species_name

    μ[node.index] = data[tip_label]
    V[node.index] = 0.0
    #V = 0.0 ## assume no measurement error
    log_nf = 0.0

    return(log_nf)
end

export simulate

function simulate(tree::Root, model::Brownian)
    data = Dict{String,Real}()

    starting_value = model.mean

    simulate_recursive(tree, model, starting_value, data)

    return(data)
end

function simulate_recursive(
        node::T, 
        model::Brownian, 
        value::Float64,
        data::Dict
    ) where {T <: InternalNode}

    left_branch = node.left
    left_node = left_branch.outbounds
    left_bl = sum(left_branch.times)
    left_d = Distributions.Normal(zero(Real), sqrt(model.sigma2 * left_bl))
    left_value = value + rand(left_d)


    right_branch = node.right
    right_node = right_branch.outbounds
    right_bl = sum(right_branch.times)
    right_d = Distributions.Normal(zero(Real), sqrt(model.sigma2 * right_bl))
    right_value = value + rand(right_d)

    
    simulate_recursive(left_node, model, left_value, data)
    simulate_recursive(right_node, model, right_value, data)
end

function simulate_recursive(
    node::Tip, 
    model::Brownian, 
    value::Float64,
    data::Dict
)
    species_name = node.species_name

    observation_distribution = Distributions.Normal(zero(Real), sqrt(model.observation_variance))
    observation_error = rand(observation_distribution)
    data[species_name] = value + observation_error
end