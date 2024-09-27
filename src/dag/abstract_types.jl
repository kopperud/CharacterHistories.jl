export DagNode

abstract type AbstractMove end
abstract type DagNode end

abstract type AbstractSet end
abstract type Continuous <: AbstractSet end
abstract type Integer <: AbstractSet end

abstract type AbstractArity end
abstract type Univariate <: AbstractArity end
abstract type Multivariate <: AbstractArity end

abstract type Stochastic <: DagNode end
abstract type Deterministic <: DagNode end
