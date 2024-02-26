export writesimmap
export simmap

function writesimmap(filename::String, tree::Root)
    s = simmap(tree)    

    open(filename, "w") do io
        write(io, newick_string)
        write(io, "\n")
    end
end

function simmap(tree::Root)
    items = String[]

    string_simmap!(tree, items)

    s = join(items)
    return(s)
end

function string_simmap!(
    node::Root,
    s::Vector{String}
)
    append!(s, ["("])

    left_node = node.left.outbounds
    right_node = node.right.outbounds

    string_simmap!(left_node, s)
    append!(s, [","])
    string_simmap!(right_node, s)

    append!(s, [");"])
end

function string_simmap!(
    node::Node,
    s::Vector{String}
)
    append!(s, ["("])

    left_node = node.left.outbounds
    right_node = node.right.outbounds

    string_simmap!(left_node, s)
    append!(s, [","])
    string_simmap!(right_node, s)

    parent_branch = node.inbounds

    append!(s, [")"])

    l = String[]
    for (state, time) in zip(parent_branch.states, parent_branch.times)
        append!(l, [join([state, time], ",")])
    end
    r = string(":{", join(l, ":"), "}")

    append!(s, [")"])

    append!(s, [r])
end


function string_simmap!(
    node::Tip,
    s::Vector{String}
)

    species_name = node.species_name  
    parent_branch = node.inbounds

    append!(s, [species_name])

    #append!(s, [")"])

    l = String[]
    for (state, time) in zip(parent_branch.states, parent_branch.times)
        append!(l, [join([state, time], ",")])
    end
    r = string(":{", join(l, ":"), "}")

    append!(s, [r])    
end