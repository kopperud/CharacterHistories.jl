export gradient_descent

## unconstrained gradient descent
function gradient_descent(
        f::Function, 
        x::Vector{Float64}; 
        max_iters::Int64 = 1000,
        ϵ::Float64 = 1e-05,
        step_size::Float64 = 1.0
    )

    objective = Inf

    state = deepcopy(x)

    global stop = false

    i = 1
    j = 0
    while !stop
        𝛁 = ForwardDiff.gradient(f, state)
        global new_state = state .- step_size .* 𝛁
        new_objective = f(new_state)

        if new_objective > objective
            step_size = step_size * 0.5
            j += 1
        else
            objective = new_objective
            state = new_state
        end

        mean_abs_gradient = sum(abs.(𝛁))/(length(𝛁))
        if mean_abs_gradient < ϵ
            println("converged in $i iters")
            stop = true
        end

        if i > max_iters
            println("did not converge in $i iters")
            stop = true
        end

        i += 1
    end
    println("halved step size $j times")
    return(objective, state)
end