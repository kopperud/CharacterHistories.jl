export Brownian

struct Brownian
    sigma2::Real
    mean::Real
    observation_variance::Real
end

Brownian(sigma2, mean) = Brownian(sigma2, mean, 0.0)

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
