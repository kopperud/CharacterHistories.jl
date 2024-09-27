export CharacterHistory, transition_probability

mutable struct CharacterHistory{T1 <: DagNode, T2 <: DagNode} <: Stochastic
    index::Int64
    value::Dict{String, String}
    tree::T1
    ch::Root
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

    CharacterHistory(dag, x1, x2, state_space)
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

#function gettree(ch::CharacterHistory)
#return(ch.tree)
#end

function parent_nodes(node::CharacterHistory)
    p = [node.tree, node.α]
    return(p)
end


function logpdf(ch::CharacterHistory)
    tree = ch.ch
   
    α = getvalue(ch.α)
    k = ch.k

    n_branches = number_of_edges(tree)
    D = zeros(eltype(α), n_branches, 2, k)

    x, log_nf = postorder!(tree, ch, D)

    ## what are appropriate root freqs?
    ## the ratio of the observed character states?
    root_freqs = [1/k for _ in 1:k]

    logl = log(sum(root_freqs .* x)) + log_nf
    return(logl)
    #return(logl, D)
end


function Base.Multimedia.display(node::CharacterHistory)
    #value = getvalue(node)
    α_value = getvalue(node.α)

    println("A Markov k-state tree distribution with rate α = $(α_value). This node has $(length(node.children)) children.")
end

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

