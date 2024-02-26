##
## in case a proposal was rejected, we need to revert back to the old state

export store
export reassign!

function store(branch::Branch)
    old_states = copy(branch.states)
    old_times = copy(branch.times)

    return(old_states, old_times)
end

function store(node::Root)
    old_left = store(node.left)
    old_right = store(node.right)

    old_state = node.state
    branches = [old_left, old_right]

    return(old_state, branches)
end

function store(node::Node)
    old_inbounds = store(node.inbounds)
    old_left = store(node.left)
    old_right = store(node.right)

    old_state = node.state
    branches = [old_inbounds, old_left, old_right]

    return(old_state, branches)
end

function reassign!(
        branch::Branch,
        old_states::Vector{String}, 
        old_times::Vector{Float64}
    )
    branch.states = old_states
    branch.times = old_times
end

function reassign!(
        node::Root,
        state, 
        branches
    )
    left, right = branches

    reassign!(node.left, left[1], left[2])
    reassign!(node.right, right[1], right[2])
    
    node.state = state
end

function reassign!(
    node::Node,
    state, 
    branches
)
    inbounds, left, right = branches

    reassign!(node.inbounds, inbounds[1], inbounds[2])
    reassign!(node.left, left[1], left[2])
    reassign!(node.right, right[1], right[2])

    node.state = state
end