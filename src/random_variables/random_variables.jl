export DagNode

abstract type DagNode end

abstract type AbstractSet end
abstract type Continuous <: AbstractSet end
abstract type Integer <: AbstractSet end

abstract type AbstractArity end
abstract type Univariate <: AbstractArity end
abstract type Multivariate <: AbstractArity end

abstract type Stochastic <: DagNode end
abstract type Deterministic <: DagNode end
#abstract type Constant <: DagNode end

export ConstantNode

mutable struct ConstantNode{T} <: DagNode
    value::T
    children::Vector{DagNode}
end

function ConstantNode(x::Float64)
    n = ConstantNode(x, DagNode[])
    return(n)
end

function ConstantNode(x::Int64)
    n = ConstantNode(x, DagNode[])
    return(n)
end

function Base.Multimedia.display(node::ConstantNode)
    value = getvalue(node)

    println("A constant node with value $(value). This node has $(length(node.children)) children.")
end


export getvalue
function getvalue(node::ConstantNode)
    return(node.value)
end



#=
abstract type UnivariateNode <: DagNode end
abstract type MultivariateNode <: DagNode end
abstract type TreeNode <: DagNode end

abstract type ContinuousUnivariateNode <: UnivariateNode end
abstract type IntegerUnivariateNode <: UnivariateNode end

abstract type RandomContinuousUnivariateNode <: ContinuousUnivariateNode end
abstract type ConstantContinuousUnivariateNode <: ContinuousUnivariateNode end


#abstract type UnivariateConstant <: ConstantVariable end
#abstract type UnivariateNatural <: ConstantVariable end



export logpdf, redraw!

export getvalue 

function getvalue(node::UnivariateConstant)
    return(node.value)
end

function UnivariateConstant(value::Float64)
    node = UnivariateConstant()
    node.value = value
    return(node)
end

#=
struct DeterministicNode
    parent::DagNode
    children::Vector{DagNode}
    x::Float64
end
=#

#update()
=#