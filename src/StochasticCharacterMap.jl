function waiting_time(model::Mk, state::String)
    d = Distributions.Exponential(1 / model.α)
    r = rand(d, 1)[1]
    return(r)
end

function sample_new_state(model::Mk, state::String)
    x = zeros(model.k)

    for (i, s) in enumerate(model.state_space)
        if s != state
            x[i] = 1/(model.k-1)
        end
    end
    d = Distributions.Categorical(x)
    idx = rand(d)
    new_state = model.state_space[idx]

    return(new_state)
end

export sample_branch!


function ancestral_character_map(tree)
    

end



function sample_branch!(
    branch::Branch, 
    model::Mk
    )
    oldest_state = branch.states[end]
    youngest_state = branch.states[1]
    states, times = sample_branch(branch, model, oldest_state, youngest_state)

    branch.states = states
    branch.times = times
end

function sample_branch!(
    branch::Branch, 
    model::Mk,
    oldest_state::String,
    youngest_state::String
    )
    states, times = sample_branch(branch, model, oldest_state, youngest_state)
    branch.states = states
    branch.times = times
end

function sample_branch(
    branch::Branch, 
    model::Mk, 
    oldest_state::String,
    youngest_state::String)

    y_state = "123868asdguu123"

    i = 1
    while youngest_state != y_state
        global states, times = sample_branch_history(branch, model, oldest_state, youngest_state)

        y_state = states[1]
        
        i += 1
        if i > 500
            throw("too many rejections, stopping at 500")
        end
    end
    return(states, times)
end

function sample_branch_history(
    branch::Branch, 
    model::Mk, 
    oldest_state::String,
    youngest_state::String)

    bl = sum(branch.times)
    states = String[]
    state_times = Float64[]

    t = [0.0]
    t_sum = 0.0
    state = oldest_state ## from old to young
    
    i = 1

    ## Nielsen (2001) Genetics, equation A2
    if oldest_state != youngest_state
        ## first waiting time
        tw = - log(1 - rand()*(1-exp(-model.α*bl)))/model.α
        if tw > bl
            throw("error")
        end
        append!(states, [state])
        append!(state_times, tw)
        append!(t, tw)
        state = sample_new_state(model, state)

        i += 1
    end

    while t_sum < bl
        tw = waiting_time(model, state)

        t_sum += tw
        if t_sum > bl
            t_sum = bl
            tw = min(bl, t_sum) - sum(t[1:i])
        end

        append!(states, [state])
        append!(state_times, tw)
        append!(t, tw)

        state = sample_new_state(model, state)
        i += 1

        if i > 50
            throw("error, too many substitutions")
        end
    end
    reverse!(states)
    reverse!(state_times)

    ## returns in the order of most recent to oldest
    return(states, state_times)
end



#function sample_branch_unequal(branch, model, oldest_state, youngest_state)