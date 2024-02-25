export loglikelihood, postorder

function loglikelihood(tree, model)
    n_branches = number_of_edges(tree)
    D = zeros(n_branches, 2, model.k)

    x, log_nf = postorder!(tree, model, D)

    ## what are appropriate root freqs?
    ## the ratio of the observed character states?
    root_freqs = [1/model.k for _ in 1:model.k]

    logl = log(sum(root_freqs .* x)) + log_nf
    return(logl, D)
end
