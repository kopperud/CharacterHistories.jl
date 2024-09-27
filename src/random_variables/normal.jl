export Normal
export parent_nodes

#mutable struct Normal <: UnivariateContinuousRV
mutable struct Normal{T1 <: DagNode, T2 <: DagNode} <: Stochastic
    index::Int64
    value::Float64 ## current state
    mu::T1 ## the mean
    sigma::T2 ## the standard deviation
    children::Vector{DagNode} ## the list of descendant nodes
end

function Normal(dag::Dag, mu::Float64, sigma::Float64)
    x1 = ConstantNode(dag, mu)
    x2 = ConstantNode(dag, sigma)
    Normal(dag, x1, x2)
end

function Normal(dag::Dag, mu::T1, sigma::T2) where {T1 <: DagNode, T2 <: DagNode}
    dag.node_counter += 1
    index = dag.node_counter

    d = Distributions.Normal(getvalue(mu), getvalue(sigma))

    value = rand(d)
    node = Normal(index, value, mu, sigma, DagNode[])
    push!(dag.nodes, node)

    push!(mu.children, node)
    push!(sigma.children, node)

    return(node)
end

export redraw!
export logpdf

function logpdf(node::Normal)
    mu_value = getvalue(node.mu)
    sigma_value = getvalue(node.sigma)
    
    d = Distributions.Normal(mu_value, sigma_value)
    lp = Distributions.logpdf(d, node.value)
    return(lp)
end

function redraw!(node::Normal)
    mu_value = getvalue(node.mu)
    sigma_value = getvalue(node.sigma)
    
    d = Distributions.Normal(mu_value, sigma_value)
    node.value = rand(d)
end

function getvalue(node::Normal)
    return(node.value)
end

function setvalue!(node::Normal, value::Float64)
    node.value = value;
end

function Base.Multimedia.display(node::Normal)
    value = getvalue(node)
    mu_value = getvalue(node.mu)
    sigma_value = getvalue(node.sigma)

    println("A Normal distribution with value $(value), mean $(mu_value) and sigma $(sigma_value). This node has $(length(node.children)) children.")
end

function parent_nodes(node::Normal)
    p = [node.mu, node.sigma]
    return(p)
end


