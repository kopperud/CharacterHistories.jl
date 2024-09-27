export Slide
export propose

mutable struct Slide <: AbstractMove
    node::Stochastic
    lambda::Float64 ## tuning parameter
    weight::Int64 
    tune_target::Float64
    acceptance_history::BitVector
end

function Slide(
        dag::Dag,
        node::T;
        lambda = 1.0, 
        weight = 1, 
        tune_target = 0.44
    ) where {T <: Stochastic}
    acceptance_history = BitVector()
    
    mv = Slide(node, lambda, weight, tune_target, acceptance_history)
    push!(dag.moves, mv)

    return(mv)
end

function tune!(move::Slide)
    ah = move.acceptance_history[end-95:end]

    @assert length(ah) > 5
    r = sum(ah)/length(ah)
    target = move.tune_target

    if r > target
        move.lambda = move.lambda*(1.0 + (r - target)/(1 - target))
    else
        move.lambda = move.lambda * target / (2*target - r)
    end

    ## delete acceptance history
    #move.acceptance_history = BitVector()

    nothing
end



function propose(dag::Dag, move::Slide) 
    node = move.node
    old_value = getvalue(node)

    old_logP = calculate_posterior(dag) 

    λ = move.lambda

    new_value = old_value + λ * (rand() -0.5)
    setvalue!(node, new_value)

    new_logP = calculate_posterior(dag)

    a1 = exp(new_logP - old_logP) 
    a2 = 1
    a = a1 * a2

    r = rand()
    accept = r < a

    if !accept
        setvalue!(node, old_value)
    end
    
    push!(move.acceptance_history, accept)

    nothing
end
