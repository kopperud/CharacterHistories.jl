#######################################
##
##               preorder
##
########################################
function preorder!(node::T, model::Mk, D::Array{Float64, 3}, F::Array{Float64, 3}, S_parent) where {T <: InternalNode}
    left_branch = node.left
    left_node = left_branch.outbounds
    left_bl = sum(left_branch.times)
    P_left = transition_probability(model, left_bl)
    F_left = S_parent ./ D[left_branch.index,2,:]
    F[left_branch.index,2,:] = F_left
    F[left_branch.index,1,:] = P_left * F_left
    S_left = F[left_branch.index,1,:] .* D[left_branch.index,1,:]
    S_left = S_left ./ sum(S_left)

    right_branch = node.right
    right_node = right_branch.outbounds
    right_bl = sum(right_branch.times)
    P_right = transition_probability(model, right_bl)
    F_right = S_parent ./ D[right_branch.index,2,:]
    F[right_branch.index,2,:] = F_right
    F[right_branch.index,1,:] = P_right * F_right
    S_right = F[right_branch.index,1,:] .* D[right_branch.index,1,:]
    S_right = S_right ./ sum(S_right)

    preorder!(left_node, model, D, F, S_left)
    preorder!(right_node, model, D, F, S_right)
end

function preorder!(node::Tip, model::Mk, D::Array{Float64, 3}, F::Array{Float64, 3}, S::Vector{Float64})
    ## do nothing
end
