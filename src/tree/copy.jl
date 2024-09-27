## shallow copy -- gives a copy of the branch, not recursively the whole tree
## don't know why I would copy this, makes a mess

function Base.copy(branch::Branch)
    branch_copy = Branch()

    branch_copy.inbounds = branch.inbounds
    branch_copy.outbounds = branch.outbounds
    branch_copy.index = branch.index
    branch_copy.states = copy(branch.states)
    branch_copy.times = copy(branch.times)

    return(branch_copy)
end

function Base.copy(node::Node)
    node_copy = Node()

    node_copy.inbounds = node.inbounds
    node_copy.left = node.left
    node_copy.right = node.right

    node_copy.index = node.index
    node_copy.state = node.state

    return(node_copy)
end

function Base.copy(node::Root)
    node_copy = Root()

    node_copy.left = node.left
    node_copy.right = node.right

    node_copy.index = node.index
    node_copy.state = node.state

    return(node_copy)
end

function Base.copy(tip::Tip)
    tip_copy = Tip()

    tip_copy.index = tip.index
    tip_copy.state = tip.state
    tip_copy.species_name = tip.species_name ## not a copy, but it's immutable anyway

    return(tip_copy)
end

export copy_node_and_connected_branches

function copy_node_and_connected_branches(node::Node)
    in = copy(node.inbounds)
    left = copy(node.left)
    right = copy(node.right)

    this_node = copy(node)

    branches = Branch[in, left, right]
    return(this_node, branches)
end

#save_branch()