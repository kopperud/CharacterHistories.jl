function hmc_iteration(
    tree, 
    data, 
    discrete_model, 
    continuous_model,
    ϵ, ## the step size 
    L, ## number of leapfrog steps
    M  ## diagonal mass matrix
    )

    M_inv = 1/M
    d = length(parameters)

    rv = Distributions.Normal(0, sqrt(M))
    phi = rand(rv)
    th_old = th

    log_p_old = NaN
end