export loglikelihood
export OrnsteinUhlenbeck

struct OrnsteinUhlenbeck{T <: Real}
    sigma2::T
    α::T
    θ::T
    observation_variance::T
end

function OrnsteinUhlenbeck(
    sigma2::T,
    α::T,
    θ::T
    ) where {T <: Real}
    x = OrnsteinUhlenbeck(sigma2, α, θ, zero(T))
    return(x)
end

function loglikelihood(
        tree::Root, 
        model::OrnsteinUhlenbeck,
        data::Dict{String, T}
    ) where {T <: Real}
    n_nodes = number_of_nodes(tree)

    μ = zeros(T, n_nodes)
    V = zeros(T, n_nodes)

    log_nf_factor = loglikelihood_po(tree, model, μ, V, data)

    ## root treatment
    ## assume that x_0 = μ_root
    root_index = tree.index
    #log_likelihood = log_nf_factor - 0.5 * log(pi * V[root_index])
    stationary_var = model.sigma2 / (2*model.α)
    d = Distributions.Normal(μ[root_index], sqrt(V[root_index] + stationary_var))
    log_likelihood = log_nf_factor + Distributions.logpdf(d, model.θ)
    return(log_likelihood)#, μ, V)
end

function loglikelihood_po(
        node::N,
        model::OrnsteinUhlenbeck, 
        μ::Vector{T},
        V::Vector{T}, 
        data::Dict{String, T}
    ) where {N <: InternalNode, T <: Real}
    α = model.α
    θ = model.θ
    σ2 = model.sigma2


    left_branch = node.left
    left_node = left_branch.outbounds
    left_bl = sum(left_branch.times)

    right_branch = node.right
    right_node = right_branch.outbounds
    right_bl = sum(right_branch.times)


    left_index = get_index(left_node)
    right_index = get_index(right_node)

    log_nf_left = loglikelihood_po(left_node, model, μ, V, data)
    log_nf_right = loglikelihood_po(right_node, model, μ, V, data)

    μ_left = μ[left_index]
    V_left = V[left_index]

    μ_left = exp(α*left_bl) * (μ_left - θ) + θ
    V_left = V_left*exp(2*α*left_bl) + (σ2/(2*α))*expm1(2*α*left_bl)

    μ_right = μ[right_index]
    V_right = V[right_index]

    μ_right = exp(α*right_bl) * (μ_right - θ) + θ
    V_right = V_right*exp(2*α*right_bl) + (σ2/(2*α))*expm1(2*α*right_bl)

    log_nf_left += α * left_bl
    log_nf_left += α * right_bl

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
    log_nf += log_nf_left #+ model.α * left_bl
    log_nf += log_nf_right #+ model.α * right_bl

    return(log_nf)
end

function loglikelihood_po(
    node::Tip, 
    model::OrnsteinUhlenbeck, 
    μ::Vector{T},
    V::Vector{T}, 
    data::Dict{String, T}
    ) where {T <: Real}
    tip_label = node.species_name

    μ[node.index] = data[tip_label]
    V[node.index] = model.observation_variance
    #V = 0.0 ## assume no measurement error
    log_nf = 0.0

    return(log_nf)
end

export loglikelihood_vcv

function loglikelihood_vcv(
        tree::Root, 
        model::OrnsteinUhlenbeck,
        data::Dict{String, T}
    ) where {T <: Real}

    tips = tip_nodes(tree)
    n_tips = number_of_tips(tree)
    nd = node_depths(tree)

    V = zeros(n_tips, n_tips)
    for i in 1:n_tips, j in 1:n_tips
        mrca1 = mrca_idx(tree, tips[j], tips[i])
        ta = nd[mrca1] ## time from root to mrca
        tij = nd[j] + nd[i] - 2*ta

        ## if root is fixed
        println("fixed root")
        #V[i,j] = (model.sigma2 / (2 * model.α)) * exp(-model.α*tij) * (1 - exp(-2*model.α*ta))
        ## If root is randomly distribued
        V[i,j] = (model.sigma2 / (2* model.α)) * exp(-model.α * tij)
    end

    #return(V)

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

function simulate(tree::Root, model::OrnsteinUhlenbeck)
    data = Dict{String,Real}()

    starting_value = model.θ

    simulate_recursive(tree, model, starting_value, data)

    return(data)
end

function simulate_recursive(
        node::T, 
        model::OrnsteinUhlenbeck, 
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
    model::OrnsteinUhlenbeck, 
    value::Float64,
    data::Dict
)
    species_name = node.species_name

    observation_distribution = Distributions.Normal(zero(Real), sqrt(model.observation_variance))
    observation_error = rand(observation_distribution)
    data[species_name] = value + observation_error
end