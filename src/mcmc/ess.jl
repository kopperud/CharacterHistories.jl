export pearson_correlation

function pearson_correlation(x::Vector{Float64}, y::Vector{Float64})
    n = length(x)
    @assert length(y) == n

    nom = n*sum(x[i]*y[i] for i in 1:n) - sum(x) * sum(y)

    denom1 = sqrt(n*sum(m^2 for m in x) - sum(x)^2)
    denom2 = sqrt(n*sum(m^2 for m in y) - sum(y)^2)

    res = nom / (denom1*denom2)

    return(res)
end

export ess

function ess(x, burnin = 0.1)
    x = x[(Int64(round(length(x)*burnin))+1):end]

    N = length(x)

    maxlag = minimum([500, N-2])

    #r = repeat([0], maxlag+1)
    r = zeros(maxlag+1)

    for k in 0:maxlag
        r[k+1] = pearson_correlation(x[1:(N-k)], x[(1+k):N])
    end 

    G = r[2:(maxlag+1)] + r[1:maxlag]

    tauf = -r[1] + 2 * G[1]

    for M in 1:(maxlag-1)
        if ((G[M+1] < G[M]) & (G[M+1] > 0))
            tauf = tauf + 2*G[M+1]
        else
            break
        end
    end

    tauf = maximum([tauf, 1])

    N_ess = N/tauf

    return(N_ess)
end