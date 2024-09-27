export LogNormal
export parent_nodes

mutable struct LogNormal{T1 <: DagNode, T2 <: DagNode} <: Stochastic
    index::Int64
    value::Float64 ## current state
    mu::T1 ## the log-mean
    sigma::T2 ## the sigma parameter
    children::Vector{DagNode} ## the list of descendant nodes
end

function LogNormal(dag::Dag, mu::Float64, sigma::Float64)
    x1 = ConstantNode(dag, mu)
    x2 = ConstantNode(dag, sigma)
    LogNormal(dag, x1, x2)
end

function LogNormal(dag::Dag, mu::T1, sigma::T2) where {T1 <: DagNode, T2 <: DagNode}
    dag.node_counter += 1
    index = dag.node_counter

    d = Distributions.LogNormal(getvalue(mu), getvalue(sigma))

    value = rand(d)
    node = LogNormal(index, value, mu, sigma, DagNode[])
    push!(dag.nodes, node)

    push!(mu.children, node)
    push!(sigma.children, node)

    return(node)
end

export redraw!
export logpdf

function logpdf(node::LogNormal)
    mu_value = getvalue(node.mu)
    sigma_value = getvalue(node.sigma)
    
    d = Distributions.LogNormal(mu_value, sigma_value)
    lp = Distributions.logpdf(d, node.value)
    return(lp)
end

function redraw!(node::LogNormal)
    mu_value = getvalue(node.mu)
    sigma_value = getvalue(node.sigma)
    
    d = Distributions.LogNormal(mu_value, sigma_value)
    node.value = rand(d)
end

function getvalue(node::LogNormal)
    return(node.value)
end

function setvalue!(node::LogNormal, value::Float64)
    node.value = value;
end

function Base.Multimedia.display(node::LogNormal)
    value = getvalue(node)
    mu_value = getvalue(node.mu)
    sigma_value = getvalue(node.sigma)

    println("A LogNormal distribution with value $(value), log-mean $(mu_value) and sigma $(sigma_value). This node has $(length(node.children)) children.")
end

function parent_nodes(node::LogNormal)
    p = [node.mu, node.sigma]
    return(p)
end


