
#######################################
##
##   ancestral state probabilities
##
########################################
export ancestral_state_probabilities

## these are for edge indices
function ancestral_state_probabilities(tree, model)
    logl, D = loglikelihood(tree, model)

    n_branches = number_of_edges(tree)
    F = zeros(n_branches, 2, model.k) ## 2 because beginning and end of branch

    left_branch_idx = tree.left.index
    right_branch_idx = tree.right.index

    F_root = ones(model.k)
    D_root = D[left_branch_idx,2,:] .* D[right_branch_idx,2,:]
    S_root = F_root .* D_root
    S_root = S_root ./ sum(S_root)

    preorder!(tree, model, D, F, S_root)
    
    S_branches = D[:,1,:] .* F[:,1,:]

    ####################################
    ##
    ##   convert to node indices
    ##
    ##########################
    n_nodes = number_of_nodes(tree)
    S_nodes = zeros(n_nodes, model.k)

    S_nodes[tree.index,:] .= S_root

    asp_po(tree, S_branches, S_nodes)

    return(S_nodes)
end

function asp_po(
        node::T, 
        S_branches::Matrix{Float64}, 
        S_nodes::Matrix{Float64}
    ) where {T <: InternalNode}

    left_branch_idx = node.left.index
    left_node_idx = node.left.outbounds.index
    S_nodes[left_node_idx,:] = S_branches[left_branch_idx,:]

    right_branch_idx = node.right.index
    right_node_idx = node.right.outbounds.index
    S_nodes[right_node_idx,:] = S_branches[right_branch_idx,:]

    asp_po(node.left.outbounds, S_branches, S_nodes)
    asp_po(node.right.outbounds, S_branches, S_nodes)

end

function asp_po(
        node::Tip, 
        S_branches::Matrix{Float64}, 
        S_nodes::Matrix{Float64}
    )
    node_idx = node.index
    S_nodes[node_idx,:] = S_branches[node.inbounds.index,:]
end
