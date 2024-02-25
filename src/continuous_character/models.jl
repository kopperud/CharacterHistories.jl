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
    sigma2::Vector{T}
    state_space::Vector{String}
    mean::T
    observation_variance::T
    k::Int64
end

export BrownianSD

function BrownianSD(
    sigma2::Vector{T}, 
    state_space::Vector{String}, 
    mean::T
    ) where {T <: Real}
     BrownianSD(sigma2, state_space, mean, zero(T), length(sigma2))
end
