export Normal

mutable struct Normal <: RandomVariable
    x::Float64 ## current state
    μ::Float64 ## mean
    σ::Float64 ## standard deviation
end

function Normal(μ::Float64, σ::Float64)
    d = Distributions.Normal(μ, σ)

    x = rand(d)
    return Normal(x, μ, σ)
end

function logpdf(rv::Normal)
    d = Distributions.Normal(rv.μ, rv.σ)
    lnl = Distributions.logpdf(d, rv.x)
    return(lnl)
end

function redraw!(rv::Normal)
    d = Distributions.Normal(rv.μ, rv.σ)

    rv.x = rand(d)
end