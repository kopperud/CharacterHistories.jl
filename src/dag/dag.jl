export Dag

mutable struct Dag
    nodes::Vector{DagNode}
    moves::Vector{AbstractMove}
    node_counter::Int64
end

function Dag()
    nodes = DagNode[]
    moves = AbstractMove[]

    node_counter = 0
    dag = Dag(nodes, moves, node_counter)
    return(dag)
end

export topological_sort

function topological_sort(dag)
    L = DagNode[]
    S = DagNode[node for node in dag.nodes if node isa ConstantNode]

    while !isempty(S) 
        node = pop!(S)
        push!(L, node)

        for child in node.children
            #parent_nodes = [n for n in parents(child)]
            pn = parent_nodes(child)

            if all([n in L for n in pn])
                push!(S, child)
            end
        end
    end

    order = [node.index for node in L]

    if length(L) != length(dag.nodes)
        error("something went wrong")
    end

    return(order)
end

export calculate_posterior

function calculate_posterior(dag::Dag) 
    node_order = topological_sort(dag)
    calculate_posterior(dag, node_order)
end

function calculate_posterior(dag::Dag, node_order::Vector{Int64})

    lnP = zeros(Float64, 1)

    for node_index in node_order 
        node = dag.nodes[node_index]

        if node isa Stochastic
            lnP[1] += logpdf(node)
        end
    end

    return(lnP[1])
end