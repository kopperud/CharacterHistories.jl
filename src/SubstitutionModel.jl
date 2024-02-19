mutable struct Mk
    α::Float64
    k::Int64
end

export Mk, transition_probability

function transition_probability(model::Mk, t::Float64)
    k = model.k
    α = model.α

    pii = (1/k) + ((k-1)/k)*exp(-k * α * t)
    pij = (1/k) - (1/k)*exp(-k * α * t)

    Q = zeros(k, k)
    Q[:,:] .= pij
    for i in 1:k
        Q[i,i] = pii
    end
    return(Q)
end

