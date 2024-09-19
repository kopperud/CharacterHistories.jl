export Exponential, setvalue!

#mutable struct Exponential <: UnivariateContinuousRV
mutable struct Exponential <: Stochastic
    value::Float64 ## current state
    rate::DagNode ## the rate parameter 
    children::Vector{DagNode}
end

function Exponential(rate::DagNode)
    d = Distributions.Exponential(getvalue(rate))

    value = rand(d)
    rv = Exponential(value, rate, DagNode[])
    push!(rate.children, rv)

    return(rv)
end

function logpdf(rv::Exponential)
    rate_value = getvalue(rv.rate)

    d = Distributions.Exponential(rate_value)
    lp = Distributions.logpdf(d, rv.value)

    return(lp)
end

function redraw!(node::Exponential)
    rate_value = getvalue(node.rate)
    d = Distributions.Exponential(rate_value)
    node.value = rand(d)
end


function setvalue!(rv::Normal, value::Float64)
    rv.value = value;
end

function getvalue(node::Exponential)
    return(node.value)
end

function Base.Multimedia.display(node::Exponential)
    value = getvalue(node)
    rate_value = getvalue(node.rate)

    println("An Exponential distribution with value $(value) and rate $(rate_value). This node has $(length(node.children)) children.")
end