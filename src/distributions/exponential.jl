export Exponential, setvalue!

#mutable struct Exponential <: UnivariateContinuousRV
mutable struct Exponential{T <: DagNode} <: Stochastic
    index::Int64
    value::Float64 ## current state
    rate::T ## the rate parameter 
    children::Vector{DagNode}
end

function Exponential(dag, rate::Float64)
    x1 = ConstantNode(dag, rate)
    Exponential(dag, x1)
end

function Exponential(dag::Dag, rate::DagNode)
    dag.node_counter += 1
    index = dag.node_counter

    d = Distributions.Exponential(getvalue(rate))

    value = rand(d)
    node = Exponential(index, value, rate, DagNode[])
    push!(rate.children, node)
    push!(dag.nodes, node)

    return(node)
end

function logpdf(node::Exponential)
    rate_value = getvalue(node.rate)

    d = Distributions.Exponential(rate_value)
    lp = Distributions.logpdf(d, node.value)

    return(lp)
end

function redraw!(node::Exponential)
    rate_value = getvalue(node.rate)
    d = Distributions.Exponential(rate_value)
    node.value = rand(d)
end


function setvalue!(node::Exponential, value::Float64)
    node.value = value;
end

function getvalue(node::Exponential)
    return(node.value)
end

function Base.Multimedia.display(node::Exponential)
    value = getvalue(node)
    rate_value = getvalue(node.rate)

    println("An Exponential distribution with value $(value) and rate $(rate_value). This node has $(length(node.children)) children.")
end


function parent_nodes(node::Exponential)
    p = [node.rate]
    return(p)
end