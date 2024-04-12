abstract type RandomVariable end

abstract type UnivariateRV <: RandomVariable end
#abstract type MultivariateRV <: RandomVariable end

abstract type NaturalUnivariateRV <: UnivariateRV end
abstract type RealUnivariateRV <: UnivariateRV end

export logpdf, redraw!

abstract type AbstractDagNode end

abstract type DagNode end

#struct DagNode{T1 <: AbstractDagNode, T2 <: RandomVariable} <: AbstractDagNode
#    parents::Vector{T1}
#    children::Vector{T1}
#    rv::T2
#end

#function update(node::DagNode)#
#
#end



struct ConstantNode <: RealUnivariateRV
    x::Float64
    children::Vector{DagNode}
    ConstantNode() = new()
end

function ConstantNode(x::Float64)
    node = ConstantNode()
    node.x = x
    return(node)
end

struct DeterministicNode
    parent::DagNode
    children::Vector{DagNode}
    x::Float64
end

struct StochasticNode
    parents::Vector{DagNode}
    children::Vector{DagNode}
    rv::RandomVariable
end

#update()
