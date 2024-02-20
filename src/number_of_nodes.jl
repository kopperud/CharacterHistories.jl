export number_of_nodes

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





