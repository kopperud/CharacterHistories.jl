export shitty_mcmc

function shitty_mcmc(
        tree1::Root, 
        cont_model, 
        discrete_model, 
        data; 
        n_iters = 5_000
    )
    tree = deepcopy(tree1)

    ntip = number_of_tips(tree)
    ninternal = number_of_internal_nodes(tree)
    logls = zeros(n_iters)
    αs = zeros(n_iters)
    means = zeros(n_iters)

    prior_mean = Distributions.Normal(8.0, 1.0)
    α_mean = Distributions.Exponential(0.3)
    

    l_discrete = loglikelihood(tree, discrete_model)[1]
    l = loglikelihood(tree, cont_model, data)

    prog = ProgressMeter.Progress(n_iters; desc = "Progress: ")

    i = 1
    while i <= n_iters
    #ProgressMeter.@showprogress for i in 1:n_iters
        ## transition rate proposal
        r = rand(Distributions.Uniform(-1,1))
        scale = exp(0.5*r)
        α_proposal = discrete_model.α * scale

        discrete_model_proposal = Mk(discrete_model.state_space, α_proposal)
        l_discrete_p = loglikelihood(tree, discrete_model_proposal)[1]

        a1 = exp(l_discrete_p - l_discrete)
        a2 = 1.0 ## proposal density, is it symmetric ? don't think so
        a = a1 * a2
        
        r1 = rand()

        if r1 < a
            ## accept, replace model
            l_discrete = l_discrete_p
            discrete_model = discrete_model_proposal
        else
            ## reject, do nothing
        end

        ## brownian rate proposal
        r = rand(Distributions.Uniform(-1,1))
        scale = exp(0.5*r)
        mean_proposal = cont_model.mean * scale

        cont_model_proposal = BrownianSD(cont_model.sigma2, cont_model.state_space, mean_proposal)
        lp = loglikelihood(tree, cont_model_proposal, data)

        a1 = exp(lp - l)
        a2 = 1.0 ## proposal density, is it symmetric ? don't think so
        a = a1 * a2
        
        r1 = rand()

        if r1 < a
            ## accept, replace model
            l = lp
            cont_model = cont_model_proposal
        else
            ## reject, do nothing
        end


        ## Character history proposal
        ## pick random Node
        node_index = rand((ntip+1):(ninternal+ntip))
        node = find_internal_node(tree, node_index)

        ## store node info incase rejection
        node_state, branch_states = store(node)

        ## redraw the node
        S = ancestral_state_probabilities(tree, discrete_model)
        redraw_node!(node, discrete_model, S)
        #logl_discrete = loglikelihood(tree, dicsrete)

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

        logls[i] = l_discrete + l
        αs[i] = discrete_model.α
        means[i] = cont_model.mean

        i += 1
        ProgressMeter.next!(prog)
    end
    return(tree, logls, αs, means)
end
    