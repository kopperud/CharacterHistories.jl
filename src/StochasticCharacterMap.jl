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

export ancestral_character_map

function redraw_nodes(node::T, model::Mk, S::Array{Float64, 2}) where {T <: InternalNode}
    ## redraw this node
    d = Distributions.Categorical(S[node.index,:])
    r = rand(d)

    new_state = model.state_space[r]
    node.state = new_state

    left_node = node.left.outbounds
    redraw_nodes(left_node, model, S)

    right_node = node.right.outbounds
    redraw_nodes(right_node, model, S)
end

function redraw_nodes(node::Tip, model::Mk, S::Array{Float64, 2})
    ## dont do anything
end

export stochastic_character_map 

function stochastic_character_map(tree, model)
    S = ancestral_state_probabilities(tree, model)

    tree2 = deepcopy(tree)

    redraw_nodes(tree2, model, S)

    redraw_branches(tree2, model)

    return(tree2)
end

function redraw_branches(node::Root, model::Mk)
    oldest_state = node.state

    ## left branch
    left_branch = node.left
    left_node = left_branch.outbounds
    left_state = left_node.state
    sample_branch!(left_branch, model::Mk, oldest_state, left_state)

    ## right branch
    right_branch = node.right
    right_node = right_branch.outbounds
    right_state = right_node.state
    sample_branch!(right_branch, model::Mk, oldest_state, right_state)

    redraw_branches(left_node, model)
    redraw_branches(right_node, model)
end

function redraw_branches(node::Node, model::Mk)
    this_node_state = node.state

    ## parent branch
    parent_branch = node.inbounds
    parent_node = parent_branch.inbounds
    parent_state = parent_node.state
    
    sample_branch!(parent_branch, model, parent_state, this_node_state)

    ## left branch
    left_branch = node.left
    left_node = left_branch.outbounds
    left_state = left_node.state
    sample_branch!(left_branch, model::Mk, this_node_state, left_state)

    ## right branch
    right_branch = node.right
    right_node = right_branch.outbounds
    right_state = right_node.state
    sample_branch!(right_branch, model::Mk, this_node_state, right_state)


    redraw_branches(left_node, model)
    redraw_branches(right_node, model)
end

function redraw_branches(node::Tip, model::Mk) 
    ## dont do anything
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