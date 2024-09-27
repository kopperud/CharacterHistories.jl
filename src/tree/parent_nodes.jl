export parent_nodes

function parent_nodes(node::Tip)
    x = Int64[]
    parent_node = node.inbounds.inbounds
    parent_nodes(parent_node, x)

    return(x)
end

function parent_nodes(node::Node, x::Vector{Int64})
    parent_node = node.inbounds.inbounds
    append!(x, node.index)
    parent_nodes(parent_node, x)
end

function parent_nodes(node::Root, x::Vector{Int64})
    append!(x, node.index)
end

#= export node_depths

## distance from the root until the node(s)
function node_depths(tree::Root)
    t = 0.0

    n = number_of_nodes(tree)
    depths = zeros(n)

    depths[tree.index] = 0.0

    node_depths(tree, depths, t)
    return(depths)
end

function node_depths(
        node::N, 
        depths::Vector{Float64}, 
        t::Float64
    ) where {N <: InternalNode}
    
    
    left_bl = sum(node.left.times)
    left_depth = t + left_bl
    left_node = node.left.outbounds

    depths[left_node.index] = left_depth

    right_bl = sum(node.right.times)
    right_depth = t + right_bl
    right_node = node.right.outbounds

    depths[right_node.index] = right_depth

    node_depths(left_node, depths, left_depth)
    node_depths(right_node, depths, right_depth)
end

function node_depths(
    node::Tip, 
    depths::Vector{Float64}, 
    t::Float64
)
    ## do nothing
end
 =#
export mrca
## find node that is the most recent common ancestor of tip 1 and tip 2
function mrca(
    tree::Root,
    tip1::Tip,
    tip2::Tip
    )
    if tip1.index == tip2.index
        mrca_node = tip1
    else
        p1 = parent_nodes(tip1)
        p2 = parent_nodes(tip2)

        common_nodes = intersect(p1, p2)
        
        mrca_node = find_internal_node(tree, common_nodes[1])
    end

    return(mrca_node)
end

export mrca_idx

function mrca_idx(
    tree::Root,
    tip1::Tip,
    tip2::Tip
    )
    if tip1.index == tip2.index
        idx = tip1.index
    else
        p1 = parent_nodes(tip1)
        p2 = parent_nodes(tip2)

        common_nodes = intersect(p1, p2)
        
        idx = common_nodes[1]
    end

    return(idx)
end

