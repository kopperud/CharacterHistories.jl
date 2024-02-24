export find_internal_node

function find_internal_node(node::T, idx::Int64) where {T <: InternalNode}
    left_node = node.left.outbounds
    right_node = node.right.outbounds

    if left_node.index == idx
        return(left_node)
    elseif right_node.index == idx
        return(right_node)
    else
        l = find_internal_node(left_node, idx)
        r = find_internal_node(right_node, idx)
        if l.index == idx
            return(l)
        elseif r.index == idx
            return(r)
        end
    end
    return(node)
end


function find_internal_node(node::Tip, idx::Int64)
    ## do nothing
    return(node)
end 

