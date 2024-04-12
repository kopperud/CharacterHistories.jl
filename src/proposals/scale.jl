export Slide
export tune!

mutable struct Slide <: AbstractMove
    rv::RandomVariable
    λ::Float64 ## tuning parameter
    weight::Int64
    tune_target::Float64
end

function Slide(rv::T) where {T <: RandomVariable}
    return(Slide(rv, 0.5, 1, 0.45))
end

function Slide(rv::T, weight::Int64) where {T <: RandomVariable}
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

function propose(rv::T, mv::Slide) where {T <: RandomVariable}
    #old_state = copy(x)

    
end