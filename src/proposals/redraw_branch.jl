export redraw_branch!

function redraw_branch(
    branch::Branch, 
    model::CharacterHistory, 
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

function redraw_branch!(
    branch::Branch, 
    model::CharacterHistory
    )
    oldest_state = branch.states[end]
    youngest_state = branch.states[1]
    states, times = redraw_branch(branch, model, oldest_state, youngest_state)

    branch.states = states
    branch.times = times
end

function redraw_branch!(
    branch::Branch, 
    model::CharacterHistory,
    oldest_state::String,
    youngest_state::String
    )
    states, times = redraw_branch(branch, model, oldest_state, youngest_state)
    branch.states = states
    branch.times = times
end
