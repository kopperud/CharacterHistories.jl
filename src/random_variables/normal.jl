export Normal

#mutable struct Normal <: UnivariateContinuousRV
mutable struct Normal{T1 <: DagNode, T2 <: DagNode} <: Stochastic
    value::Float64 ## current state
    mu::T1 ## the mean
    sigma::T2 ## the standard deviation
    children::Vector{DagNode} ## the list of descendant nodes
end

function Normal(mu::T1, sigma::T2) where {T1 <: DagNode, T2 <: DagNode}
    d = Distributions.Normal(getvalue(mu), getvalue(sigma))

    value = rand(d)
    rv = Normal(value, mu, sigma, DagNode[])
    push!(mu.children, rv)
    push!(sigma.children, rv)

    return(rv)
end

export redraw!
export logpdf

function logpdf(rv::Normal)
    mu_value = getvalue(rv.mu)
    sigma_value = getvalue(rv.sigma)
    
    d = Distributions.Normal(mu_value, sigma_value)
    lp = Distributions.logpdf(d, rv.value)
    return(lp)
end

function redraw!(rv::Normal)
    mu_value = getvalue(rv.mu)
    sigma_value = getvalue(rv.sigma)
    
    d = Distributions.Normal(mu_value, sigma_value)
    rv.value = rand(d)
end

function getvalue(node::Normal)
    return(node.value)
end

function setvalue!(node::Exponential, value::Float64)
    node.value = value;
end

function Base.Multimedia.display(node::Normal)
    value = getvalue(node)
    mu_value = getvalue(node.mu)
    sigma_value = getvalue(node.sigma)

    println("A Normal distribution with value $(value), mean $(mu_value) and sigma $(sigma_value). This node has $(length(node.children)) children.")
end


