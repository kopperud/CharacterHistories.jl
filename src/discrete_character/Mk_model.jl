export Mk, transition_probability

mutable struct Mk{T <: Real}
    state_space::Vector{String}
    α::T
    k::Int64
end

## easier constructor
function Mk(state_space::Vector{String}, α::T) where {T <: Real}
    Mk(state_space, α, length(state_space))
end

function transition_probability(
    model::Mk, 
    t::Float64
    )
    k = model.k
    α = model.α

    pii = (1/k) + ((k-1)/k)*exp(-k * α * t)
    pij = (1/k) - (1/k)*exp(-k * α * t)

    Q = zeros(eltype(α), k, k)
    Q[:,:] .= pij
    for i in 1:k
        Q[i,i] = pii
    end
    return(Q)
end

