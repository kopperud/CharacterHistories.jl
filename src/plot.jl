export coordinates
export treeplot

function treeplot() end

function coordinates(tree::Root)
    n_edges = number_of_edges(tree)

    y = zeros(n_edges)
    find_y!(tree.left.outbounds, y)
    find_y!(tree.right.outbounds, y)

    root_height = maximum(node_depths(tree))

    x = Dict{Int64,Vector{Float64}}()
    states = Dict{Int64, Vector{String}}()
    find_x!(tree, x, states, root_height)

    return(x, y, states)
end

function find_y!(node::Node, y::Vector{Float64})
    left_branch = node.left
    left_node = left_branch.outbounds

    right_branch = node.right
    right_node = right_branch.outbounds

    find_y!(left_node, y)
    find_y!(right_node, y)

    index = node.inbounds.index ## parental branch index 
    y[index] = (y[left_branch.index] + y[right_branch.index]) / 2
end

function find_y!(node::Tip, y::Vector{Float64})
    index = node.inbounds.index  ## parental branch index 
    
    y[index] = node.index
end

function find_x!(
    node::T,
    x::Dict{Int64,Vector{Float64}},
    states::Dict{Int64, Vector{String}},
    t::Float64
    ) where {T <: InternalNode}

    for branch in (node.left, node.right)
        bl = sum(branch.times)
        t1 = t - bl
        times = t .- cumsum(reverse(branch.times))
        prepend!(times, t)

        x[branch.index] = times
        states[branch.index] = reverse(branch.states)
        find_x!(branch.outbounds, x, states, t1)
    end
end

function find_x!(::Tip, ::Dict{Int64}, states, t::Float64)
end

export node_depths

## distance from the root to the node
function node_depths(tree::Root)
    n_nodes = number_of_nodes(tree)
    depths = zeros(Float64, n_nodes)

    t = 0.0

    node_depths_po!(tree, depths, t)    
    return(depths)
end

function node_depths_po!(node::T, depths::Vector{Float64}, t::Float64) where {T <: InternalNode}
    depths[node.index] = t

    tl = sum(node.left.times)
    node_depths_po!(node.left.outbounds, depths, t + tl)

    tr = sum(node.right.times)
    node_depths_po!(node.right.outbounds, depths, t + tr)
end

function node_depths_po!(node::Tip, depths::Vector{Float64}, t::Float64)
    depths[node.index] = t
end

