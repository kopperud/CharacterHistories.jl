export sample_branch!

function waiting_time(model::CharacterHistory, state::String)
    α = getvalue(model.α)
    d = Distributions.Exponential(1 / α)
    r = rand(d, 1)[1]
    return(r)
end

function sample_new_state_naive(model::CharacterHistory, state::String)
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

function sample_new_state(model::CharacterHistory, state::String)
    r = rand()
    p = 1.0/(model.k-1.0)

    #for (i, s) in enumerate(model.state_space)
    for i in eachindex(model.state_space)
        #if s != state
        if model.state_space[i] != state
            if r < p
                global idx = i
                break
            else
                r -= p
            end
        end
    end

    new_state = model.state_space[idx]

    return(new_state)
end

function sample_branch_history(
    branch::Branch, 
    model::CharacterHistory, 
    oldest_state::String,
    youngest_state::String)

    bl = sum(branch.times)
    states = String[]
    state_times = Float64[]

    t = [0.0]
    t_sum = 0.0
    state = oldest_state ## from old to young
    
    i = 1

    α = getvalue(model.α)

    ## Nielsen (2001) Genetics, equation A2
    if oldest_state != youngest_state
        ## first waiting time
        tw = - log(1 - rand()*(1-exp(-α*bl)))/α
        if tw > bl
            throw("error")
        end
        append!(states, [state])
        append!(state_times, tw)
        append!(t, tw)
        state = sample_new_state(model, state)
        t_sum += tw

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

        if i > 500
            throw("error, too many substitutions (rate = $(model.α)")
        end
    end
    reverse!(states)
    reverse!(state_times)

    ## returns in the order of most recent to oldest
    return(states, state_times)
end