export loglikelihood
export OrnsteinUhlenbeckSD

struct OrnsteinUhlenbeckSD{T <: Real}
    parameters::OrderedCollections.LittleDict{String, Vector{T}}
    observation_variance::T
    k::Int64
end

function OrnsteinUhlenbeckSD(
    state_space::Vector{String},
    sigma2s::Vector{T},
    αs::Vector{T},
    θs::Vector{T}
    ) where {T <: Real}
    @assert length(sigma2s) == length(αs)
    @assert length(αs) == length(θs)

    d = OrderedCollections.LittleDict{String,Vector{T}}()
    for (state, σ2, α, θ) in zip(state_space, sigma2s, αs, θs)
        d[state] = [σ2, α, θ]
    end

    OrnsteinUhlenbeckSD(d, zero(T), length(sigma2s))
end

function loglikelihood(
        tree::Root, 
        model::OrnsteinUhlenbeckSD,
        data::Dict{String, T}
    ) where {T <: Real}
    n_nodes = number_of_nodes(tree)

    μ = zeros(T, n_nodes)
    V = zeros(T, n_nodes)

    log_nf_factor = loglikelihood_po(tree, model, μ, V, data)

    ## root treatment
    ## assume that x_0 = μ_root
    root_index = tree.index
    root_state = tree.state

    root_theta = model.parameters[root_state][3]
    #log_likelihood = log_nf_factor - 0.5 * log(pi * V[root_index])
    d = Distributions.Normal(μ[root_index], sqrt(V[root_index]))
    log_likelihood = log_nf_factor + Distributions.logpdf(d, root_theta)
    return(log_likelihood)#, μ, V)
end

function loglikelihood_po(
        node::N,
        model::OrnsteinUhlenbeckSD, 
        μ::Vector{T},
        V::Vector{T}, 
        data::Dict{String, T}
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
        σ2, α, θ = model.parameters[state]
        vy = σ2 / (2*α)

        V_left = V_left * exp(2* α * time) + vy*expm1(2*α*time)
        μ_left = exp(α * time) * (μ_left - θ) + θ

        log_nf_left += time * α
    end

    right_index = get_index(right_node)

    μ_right = μ[right_index]
    V_right = V[right_index]
    for (time, state) in zip(right_branch.times, right_branch.states)
        σ2, α, θ = model.parameters[state]
        vy = σ2 / (2*α)

        V_right = V_right * exp(2 * α * time) + vy*expm1(2*α*time)
        μ_right = exp(α * time) * (μ_right - θ) + θ

        log_nf_right += time * α
    end

    ## merging rule
    μ_node = (μ_left*V_right + μ_right*V_left) / (V_left+V_right)
    V_node = (V_left*V_right)/(V_left+V_right)

    #println(μ_left)
    node_index = get_index(node)
    μ[node_index] = μ_node
    V[node_index] = V_node

    ## normalizing factor
    contrast = μ_left - μ_right
    log_nf = -(1.0)*contrast^2 / (2*(V_left+V_right))
    log_nf -= 0.5 * log(2 * pi * (V_left+V_right))
    log_nf += log_nf_left
    log_nf += log_nf_right

    return(log_nf)
end

function loglikelihood_po(
    node::Tip, 
    model::OrnsteinUhlenbeckSD, 
    μ::Vector{T},
    V::Vector{T}, 
    data::Dict{String, T}
    ) where {T <: Real}
    tip_label = node.species_name

    index = get_index(node)
    μ[index] = data[tip_label]
    V[index] = model.observation_variance
    #V = 0.0 ## assume no measurement error
    log_nf = 0.0

    return(log_nf)
end

export loglikelihood_vcv

function loglikelihood_vcv(
        tree::Root, 
        model::OrnsteinUhlenbeckSD,
        data::Dict{String, T}
    ) where {T <: Real}

    tips = tip_nodes(tree)
    n_tips = number_of_tips(tree)
    nd = node_depths(tree)

    V = zeros(n_tips, n_tips)
    for i in 1:n_tips, j in 1:n_tips
        mrca1 = mrca_idx(tree, tips[j], tips[i])
        V[i,j] = nd[mrca1] * model.sigma2
    end

    Ve = LinearAlgebra.diagm(ones(n_tips) .* model.observation_variance)
    V = V + Ve

    
    X = ones(n_tips, 1)
    tl = tiplabels(tree)
    y = [data[tl[i]] for i in 1:n_tips]
    β = [model.θ]

    logl = gls_logpdf(V, X, y, β)
    return(logl)
end

export simulate

function simulate(tree::Root, model::OrnsteinUhlenbeckSD)
    data = Dict{String,Real}()

    starting_value = model.θ

    simulate_recursive(tree, model, starting_value, data)

    return(data)
end

function simulate_recursive(
        node::T, 
        model::OrnsteinUhlenbeckSD, 
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
    model::OrnsteinUhlenbeckSD, 
    value::Float64,
    data::Dict
)
    species_name = node.species_name

    observation_distribution = Distributions.Normal(zero(Real), sqrt(model.observation_variance))
    observation_error = rand(observation_distribution)
    data[species_name] = value + observation_error
end