export Scale
export tune!

mutable struct Scale <: AbstractMove
    node::Stochastic
    lambda::Float64 ## tuning parameter
    weight::Int64 
    tune_target::Float64
    acceptance_history::BitVector
end

function Scale(
    dag::Dag,
    node::T1;
    lambda::Float64 = 1.0,
    weight::Int64 = 1,
    tune_target = 0.44,
) where {T1 <: Stochastic}
    acceptance_history = BitVector()

    mv = Scale(node, lambda, weight, tune_target, acceptance_history)
    push!(dag.moves, mv)

    return(mv)
end

function Base.Multimedia.display(move::Scale)
    lambda = move.lambda
    w = move.weight

    println("A scale move with lambda $(lambda) and weight $w.") 
end


function tune!(move::Scale)
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


function propose(dag::Dag, move::Scale)
    node = move.node
    old_value = getvalue(node)

    old_logP = calculate_posterior(dag) 

    lambda = move.lambda

    u = rand()
    sf = exp(lambda*(u-0.5))

    new_value = old_value * sf
    setvalue!(node, new_value)

    new_logP = calculate_posterior(dag)

    a1 = exp(new_logP - old_logP) 
    a2 = sf

    a = a1 * a2

    r = rand()
    accept = r < a
    if !accept
        setvalue!(node, old_value)
    end

    push!(move.acceptance_history, accept)
    
    return(nothing)
end
