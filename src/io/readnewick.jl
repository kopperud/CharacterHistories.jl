## read newick
export readfile
export tokenize
export readsimmap
export partition

function readfile(filename)
    io = open(filename, "r")
    s = read(io, String)
    close(io)

    return(s)
end

function readsimmap(filename)
    s = readfile(filename)

    s = replace(s, r"^_+|_+$" => "")
    s = replace(s, r"[ \t]" => "") ## delete comments
    s = replace(s, '\n' => "")
    #'\n' => ""
    s = s[findfirst('(', s):end]
    s = s[1:findfirst(';', s)]
    s = stripcomments(s)
    tokens = tokenize(s)

    tokens2 = tokens[2:end-2]
    left, right = partition(tokens2)

    ntip = length(findall(')', s))+1

    if left[end][1] != ':'
        left_branch = terminal(left[1])
    else
        left_branch = internal(left)
    end

    if right[end][1] != ':'
        right_branch = terminal(right[1])
    else
        right_branch = internal(right)
    end

    root = Root()
    root.left = left_branch
    root.right = right_branch

    right_branch.inbounds = root
    left_branch.inbounds = root
    root.state = root.left.states[end]

    reindex!(root) ## make indices on the tree nodes and branches

    return(root)
end

function stripcomments(s)
    res = replace(s, r"\[.*?\]" => "")
    return(res)
end

function tokenize(s)
    tokens = String[]

    ## strip everything between square bracket
    s = stripcomments(s)
    len = length(s)

    single_tokens = Set([')', '(', ',', ';']) # , '{', '}'
    i = 1
    while i <= len
        if s[i] ∈ single_tokens
            token = string(s[i])
            append!(tokens, [token])
            i += 1
        else
            ## if is a simmap comment
            if occursin(":{", s[i:end])
                close_idx = findfirst('}', @view s[i:end])
                token = @view s[i:close_idx+i-1]
                append!(tokens, [token])
                i += length(token)
            else
                l = Int64[]
                firstcomma = findfirst(',', @view s[i:end])
                firstclose = findfirst(')', @view s[i:end])
                if !isnothing(firstcomma)
                    append!(l, firstcomma-1)
                end
                if !isnothing(firstclose)
                    append!(l, firstclose-1)
                end

                if !isempty(l)
                    close_idx = minimum(l)
                    token = @view s[i:close_idx+i-1]
                    append!(tokens, [token])
                    i += length(token)
                else
                    i += 1
                end
            end
        end
    end
    return(tokens)
end

export parse_history

function parse_history(s)
    curly_idx = findfirst('{', s)
    h = s[(curly_idx+1):(end-1)]
    states = String[]
    times = Float64[]

    for item in split(h, ':')
        state, time = split(item, ',')
        append!(states, [state])
        append!(times, parse(Float64, time))
    end
    return(states, times)
end

function parse_brlen(s)
    res = parse(Float64, split(s, ':')[end])
    return(res)
end

function parse_tiplab(s)
    idx = findfirst(':', s)-1
    res = s[1:idx]
    return(res)
end

function findsplit(tokens)
    global ps = 0
    for (i, token) in enumerate(tokens)
        if token == "("
            ps += 1
        elseif token == ")"
            ps -= 1
        end
        
        if (token == ",") & (ps == 0)
            return(i)
        end
    end
    throw("split not found") 
end


function partition(tokens)
    comma = findsplit(tokens)

    left = @view tokens[1:comma-1]
    right = @view tokens[1+comma:end]

    return (left, right)
end


function internal(tokens)
    states, times = parse_history(tokens[end])
    tokens = tokens[2:end-2]

    left, right = partition(tokens)

    node = Node()

    if left[end][1] != ':'
        left_branch = terminal(left[1])
    else
        left_branch = internal(left)
    end

    if right[end][1] != ':'
        right_branch = terminal(right[1])
    else
        right_branch = internal(right)
    end
    
    node.left = left_branch
    node.right = right_branch

    right_branch.inbounds = node
    left_branch.inbounds = node

    branch = Branch()
    branch.outbounds = node
    branch.states = states
    branch.times = times

    node.inbounds = branch
    node.state = states[end]
    
    return(branch)
end

function terminal(s)
    tiplab = parse_tiplab(s)

    tip = Tip()
    tip.species_name = tiplab
    
    states, times = parse_history(s)

    branch = Branch()
    branch.outbounds = tip
    branch.states = states
    branch.times = times

    tip.inbounds = branch
    tip.state = states[1]
    
    return(branch)
    return(tip)
end