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