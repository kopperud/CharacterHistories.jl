export Slide
export tune!

mutable struct Slide <: AbstractMove
    node::Stochastic
    lambda::Float64 ## tuning parameter
    weight::Int64 
    tune_target::Float64
end

function Slide(node::T; lambda = 1.0, weight = 1, tune_target = 0.23) where {T <: Stochastic}
    return(Slide(node, lambda, weight, tune_target))
end

export propose

function propose(dag::Dag, node::T, move::Slide) where {T <: Stochastic}
    old_value = getvalue(node)

    old_logP = calculate_posterior(dag) 

    λ = move.lambda

    new_value = old_value + λ * (rand() -0.5)
    setvalue!(node, new_value)

    new_logP = calculate_posterior(dag)

    a1 = exp(new_logP - old_logP) 
    a2 = 1
    a = a1 * a2

    r = rand()
    if r > a
        setvalue!(node, old_value)
    end

    nothing
end

export Scale
mutable struct Scale 
    λ::Float64
end

function propose(dag::Dag, node::T, move::Scale) where {T <: Stochastic}
    old_value = getvalue(node)

    old_logP = calculate_posterior(dag) 

    λ = move.λ

    u = rand()
    sf = exp(λ*(u-0.5))

    new_value = old_value * sf
    setvalue!(node, new_value)

    new_logP = calculate_posterior(dag)

    a1 = exp(new_logP - old_logP) 
    a2 = sf

    a = a1 * a2

    r = rand()
    if r > a
        setvalue!(node, old_value)
    end
    
    nothing
end



#=
mutable struct Slide <: AbstractMove
    rv::UnivariateContinuousRV
    λ::Float64 ## tuning parameter
    weight::Int64
    tune_target::Float64
end

function Slide(rv::T) where {T <: UnivariateContinuousRV}
    return(Slide(rv, 0.5, 1, 0.45))
end

function Slide(rv::T, weight::Int64) where {T <: UnivariateContinuousRV}
    return(Slide(rv, 0.5, weight, 0.45))
end

function tune!(mv::Slide, acceptance::Vector{Bool}, window::Int64)
    acceptence_target = mv.tune_target
    acceptance_rate = sum(acceptance[end-window+1:end])/window

    λ = mv.λ

    if λ > acceptence_target
        λ_new = λ + λ*((acceptance_rate-acceptence_target)/(1.0-acceptence_target))
    else
        λ_new = λ / (2.0 - (acceptence_target/acceptence_target))
    end
    mv.λ = λ_new
end

function propose(rv::T, mv::Slide) where {T <: UnivariateContinuousRV}
    #old_state = copy(x)

    
end

=#
