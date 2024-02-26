export shitty_mcmc

function shitty_mcmc(
        tree1::Root, 
        cont_model, 
        discrete_model, 
        data; 
        n_iters = 5_000
    )
    tree = deepcopy(tree1)

    S = ancestral_state_probabilities(tree, discrete_model)
    ntip = number_of_tips(tree)
    ninternal = number_of_internal_nodes(tree)
    logls = zeros(n_iters)
    l = loglikelihood(tree, cont_model, data)

    ProgressMeter.@showprogress for i in 1:n_iters
        ## pick random Node
        node_index = rand((ntip+1):(ninternal+ntip))
        node = find_internal_node(tree, node_index)

        ## store node info incase rejection
        node_state, branch_states = store(node)

        ## redraw the node
        redraw_node!(node, discrete_model, S)

        lp = loglikelihood(tree, cont_model, data)
        a1 = exp(lp - l)
        a2 = 1.0 ## proposal density
        a = a1 * a2

        r1 = rand()

        if r1 < a
            ## accept, do nothing
            l = lp
        else
            ## reject, replace old values
            reassign!(node, node_state, branch_states)
        end

        logls[i] = l
    end
    return(tree, logls)
end
    