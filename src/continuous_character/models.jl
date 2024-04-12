export Brownian

struct Brownian{T <: Real}
    sigma2::T
    mean::T
    observation_variance::T
end

function Brownian(
    sigma2::T,
    mean::T
    ) where {T <: Real}
    x = Brownian(sigma2, mean, zero(T))
    return(x)
end

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


struct OrnsteinUhlenbeckSD{T <: Real}
    parameters::OrderedCollections.LittleDict{String, Vector{T}}
    observation_variance::T
    k::Int64
end

export OrnsteinUhlenbeckSD

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
