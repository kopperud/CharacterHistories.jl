export CharacterHistory, transition_probability

mutable struct CharacterHistory{T1 <: DagNode, T2 <: DagNode} <: Stochastic
    index::Int64
    value::Dict{String, String}
    tree::T1  ## this is the parent in the DAG
    ch::Root  ## this is the tree with the augmented character history
    state_space::Vector{String}
    α::T2
    k::Int64
    children::Vector{DagNode}
end

sample(rng::Random.AbstractRNG, a::Random.AbstractArray) = a[rand(rng, 1:length(a))]
sample(a::AbstractArray) = sample(Random.default_rng(), a)

function CharacterHistory(
        dag::Dag,
        tree::Root,
        α::Float64,
        state_space::Vector{String}
    )
    x1 = ConstantNode(dag, tree)
    x2 = ConstantNode(dag, α)

    CharacterHistory(dag, tree, x2, state_space)
end

## easier constructor
function CharacterHistory(
        dag::Dag, 
        tree::T1, 
        α::T2, 
        state_space::Vector{String}
    ) where {T1 <: DagNode, T2 <: DagNode}
    dag.node_counter += 1 
    index = dag.node_counter

    ch = copy(getvalue(tree))

    #tl = tiplabels(tr)
    #state_space = unique(values(tipstates(tr))) ## problem: this can reorder the tip states
    #tip_states = Dict{String,String}(species => sample(state_space) for species in tl) 

    ## todo: the input tree should not be a character history,
    ## and then this needs to be changed somehow
    value = tipstates(ch)
    
    children = DagNode[]
    node = CharacterHistory(index, value, tree, ch, state_space, α, length(state_space), children)
    push!(dag.nodes, node)

    push!(tree.children, node)
    push!(α.children, node)

    return node
end

function getvalue(ch::CharacterHistory)
    return(ch.value)
end

function get_tree(ch::CharacterHistory)
    return(ch.ch)
end

function parent_nodes(node::CharacterHistory)
    p = [node.tree, node.α]
    return(p)
end


function logpdf(ch::CharacterHistory)
    root_node = ch.ch
   
    α = getvalue(ch.α)
    k = ch.k

    n_branches = number_of_edges(root_node)
    D = zeros(eltype(α), n_branches, 2, k)

    x, log_nf = postorder!(root_node, ch, D)

    ## what are appropriate root freqs?
    ## the ratio of the observed character states?
    root_freqs = [1/k for _ in 1:k]

    logl = log(sum(root_freqs .* x)) + log_nf
    return(logl)
    #return(logl, D)
end


function Base.Multimedia.display(node::CharacterHistory)
    α_value = getvalue(node.α)

    println("A Markov k-state tree distribution with rate α = $(α_value). This node has $(length(node.children)) children.")
end

#######################################
##
##   transition probabilities, i.e. the P(Δt) matrix 
##
########################################
function transition_probability(
    ch::CharacterHistory, 
    #α::Float64,
    t::Float64
    )
    k = ch.k
    α = getvalue(ch.α)

    pii = (1/k) + ((k-1)/k)*exp(-k * α * t)
    pij = (1/k) - (1/k)*exp(-k * α * t)

    Q = zeros(eltype(α), k, k)
    Q[:,:] .= pij
    for i in 1:k
        Q[i,i] = pii
    end
    return(Q)
end


#######################################
##
##               postorder
##
########################################
function postorder!(
    node::N, 
    ch::CharacterHistory,
    D::Array{T, 3}
    ) where {N <: InternalNode, T <: Real}

    #α = getvalue(ch.α)

    left_branch = node.left
    left_node = left_branch.outbounds
    left_bl = sum(left_branch.times)
    P_left = transition_probability(ch, left_bl)

    right_branch = node.right
    right_node = right_branch.outbounds
    right_bl = sum(right_branch.times)
    P_right = transition_probability(ch, right_bl)

    x_left, log_nf_left = postorder!(left_node, ch, D)
    x_right, log_nf_right = postorder!(right_node, ch, D)

    D[left_branch.index,1,:] = x_left
    D[right_branch.index,1,:] = x_right
    
    state_left = P_left * x_left
    state_right = P_right * x_right

    D[left_branch.index,2,:] = state_left
    D[right_branch.index,2,:] = state_right
    
    state = state_left .* state_right
    nf_node = sum(state)
    state = state ./ nf_node
    log_nf = log_nf_left + log_nf_right + log(nf_node)

    return(state, log_nf)
end

function postorder!(
        node::Tip, 
        ch::CharacterHistory,
        D::Array{T, 3}
    ) where {T <: Real}

    state = zeros(ch.k)
    idx = findfirst(isequal(node.state), ch.state_space)
    state[idx] = 1
    log_nf = 0.0

    return(state, log_nf)
end

#######################################
##
##               preorder
##
########################################
function preorder!(
        node::T, 
        model::CharacterHistory, 
        D::Array{Float64, 3}, 
        F::Array{Float64, 3}, 
        S_parent::Vector{Float64}
    ) where {T <: InternalNode}
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

function preorder!(
    node::Tip, 
    model::CharacterHistory, 
    D::Array{Float64, 3}, 
    F::Array{Float64, 3}, 
    S::Vector{Float64}
    )
    ## do nothing
    nothing
end


#######################################
##
##   ancestral state probabilities
##
########################################
export ancestral_state_probabilities

## these are for edge indices
function ancestral_state_probabilities(
    ch::CharacterHistory,
)
    α = getvalue(ch.α)
    k = ch.k
    root = ch.ch
    n_branches = number_of_edges(root)
    F = zeros(eltype(α), n_branches, 2, k) ## 2 because beginning and end of branch
    D = zeros(eltype(α), n_branches, 2, k)

    ## Postorder    
    x, log_nf = postorder!(root, ch, D) ## this fills in D

    ## some shenanigans at the root
    left_branch_idx = root.left.index
    right_branch_idx = root.right.index

    F_root = ones(k)
    D_root = D[left_branch_idx,2,:] .* D[right_branch_idx,2,:]
    S_root = F_root .* D_root
    S_root = S_root ./ sum(S_root)

    ## Preorder
    preorder!(root, ch, D, F, S_root)
   
    ## ancestral state probs
    S_branches = D[:,1,:] .* F[:,1,:]

    ####################################
    ##
    ##   convert to node indices
    ##
    ##########################
    n_nodes = number_of_nodes(root)
    S_nodes = zeros(n_nodes, k)

    S_nodes[root.index,:] .= S_root
    asp_po(root, S_branches, S_nodes)

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
