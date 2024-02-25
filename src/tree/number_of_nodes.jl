export number_of_nodes

## number of nodes
function number_of_nodes(tree::Root)
    n = 0
    n = nnodes_po(tree, n)
    return(n)
end

function nnodes_po(::Tip, n::Int64)
    return(n+1)
end

function nnodes_po(node::T, n::Int64) where {T <: InternalNode}
    n = nnodes_po(node.left.outbounds, n)
    n = nnodes_po(node.right.outbounds, n)
    return(n+1)
end

## number of internal nodes
export number_of_internal_nodes
function number_of_internal_nodes(tree::Root)
    n = 0
    n = n_internal_po(tree, n)
    return(n)
end

function n_internal_po(node::T, n::Int64) where {T <: InternalNode}
    n = n_internal_po(node.left.outbounds, n)
    n = n_internal_po(node.right.outbounds, n)
    return(n+1)
end
function n_internal_po(::Tip, n::Int64) return(n) end

export number_of_edges

function number_of_edges(tree::Root)
    n = 0
    n = nedges_po(tree, n)
    return(n)
end

function nedges_po(::Tip,n::Int64) 
    return(n)
end

function nedges_po(node::T, n::Int64) where {T <: InternalNode}
    n = nedges_po(node.left.outbounds, n)
    n = nedges_po(node.right.outbounds, n)
    return(n+2)
end



export number_of_tips

function number_of_tips(tree::Root)
    i = zeros(Int64, 1)
    number_of_tips(tree, i)
    return(i[1])
end

function number_of_tips(
    node::N,
    i::Vector{Int64}
) where {N <: InternalNode}
    left_node = node.left.outbounds
    right_node = node.right.outbounds

    number_of_tips(left_node, i)
    number_of_tips(right_node, i)
end

function number_of_tips(
    node::Tip,
    i::Vector{Int64}
)
    i[1] += 1
end



