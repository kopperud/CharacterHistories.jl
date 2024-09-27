export ConstantNode

mutable struct ConstantNode{T} <: DagNode
    index::Int64
    value::T
    children::Vector{DagNode}
end

function ConstantNode(dag::Dag, value::Root)
    dag.node_counter += 1
    index = dag.node_counter

    node = ConstantNode(index, value, DagNode[])
    push!(dag.nodes, node)

    return(node)
end

function ConstantNode(dag::Dag, value::Float64)
    dag.node_counter += 1
    index = dag.node_counter

    node = ConstantNode(index, value, DagNode[])
    push!(dag.nodes, node)

    return(node)
end

function ConstantNode(dag::Dag, value::Vector{Float64})
    dag.node_counter += 1
    index = dag.node_counter

    node = ConstantNode(index, value, DagNode[])
    push!(dag.nodes, node)

    return(node)
end

function ConstantNode(dag::Dag, value::Matrix{Float64})
    dag.node_counter += 1
    index = dag.node_counter

    node = ConstantNode(index, value, DagNode[])
    push!(dag.nodes, node)

    return(node)
end



function Base.Multimedia.display(node::ConstantNode)
    value = getvalue(node)

    println("A constant node with value $(value). This node has $(length(node.children)) children.")
end

function Base.Multimedia.display(node::ConstantNode{Root})
    value = getvalue(node)

    println("A constant node representing a tree with a history of discrete character states. This node has $(length(node.children)) children.")
end



export getvalue
function getvalue(node::ConstantNode)
    return(node.value)
end

