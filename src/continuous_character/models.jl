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
