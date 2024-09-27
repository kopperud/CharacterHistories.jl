export loglikelihood

const PI = 3.14159265358979


struct BrownianSD{T <: Real}
    sigma2::OrderedCollections.LittleDict{String, T}
    #state_space::Vector{String}
    mean::T
    observation_variance::T
    k::Int64
end

export BrownianSD

function BrownianSD(
    state_space::Vector{String},
    sigma2::Vector{T},
    mean::T
    ) where {T <: Real}
    d = OrderedCollections.LittleDict{String,T}()
    for (state, s2) in zip(state_space, sigma2)
        d[state] = s2
    end

    BrownianSD(d, mean, zero(T), length(sigma2))
end



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


export simulate

function simulate(tree::Root, model::BrownianSD)

    sample_sizes = Dict{String,Int64}()
    species_names = tiplabels(tree)
    for name in species_names
        sample_sizes[name] = 1
    end

    data = Dict{String,Real}()
    starting_value = model.mean

    simulate_recursive(tree, model, starting_value, data, sample_sizes)

    return(data)
end

function simulate(tree::Root, model::BrownianSD, sample_sizes::Dict)
    data = Dict{String,Real}()
    starting_value = model.mean

    simulate_recursive(tree, model, starting_value, data, sample_sizes)

    return(data)
end

function simulate_recursive(
        node::T, 
        model::BrownianSD, 
        value::Float64,
        data::Dict,
        sample_sizes::Dict
    ) where {T <: InternalNode}

    left_branch = node.left
    left_node = left_branch.outbounds
    left_value = value

    times = reverse(left_branch.times)
    states = reverse(left_branch.states)
    for (time, state) in zip(times, states)
        left_d = Distributions.Normal(zero(Real), sqrt(model.sigma2[state] * time))
        left_value += rand(left_d)
    end      


    right_branch = node.right
    right_node = right_branch.outbounds
    right_value = value

    times = reverse(right_branch.times)
    states = reverse(right_branch.states)
    for (time, state) in zip(times, states)
        right_d = Distributions.Normal(zero(Real), sqrt(model.sigma2[state] * time))
        right_value += rand(right_d)
    end  
    
    simulate_recursive(left_node, model, left_value, data, sample_sizes)
    simulate_recursive(right_node, model, right_value, data, sample_sizes)
end

function simulate_recursive(
    node::Tip, 
    model::BrownianSD, 
    value::Float64,
    data::Dict,
    sample_sizes::Dict
)
    species_name = node.species_name
    n = sample_sizes[species_name]

    observation_distribution = Distributions.Normal(zero(Real), sqrt(model.observation_variance))
    observation_error = rand(observation_distribution)
    data[species_name] = value + observation_error
end