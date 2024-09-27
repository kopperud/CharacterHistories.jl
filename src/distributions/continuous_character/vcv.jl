export gls_logpdf

function gls_logpdf(
    V::Matrix{Float64},
    X::Matrix{Float64},
    y::Vector{Float64},
    β::Vector{Float64}
)
    C = LinearAlgebra.cholesky(V)
    L = C.L

    n = length(y)

    log_det_V = sum([log(L[i,i]) for i in 1:n])*2.0

    #y = [data[species] for species in tree[:tip_label]]
    r = L\y .- L\(X * β) ## de-correlated residuals

    res = 0.0
    res -= (n/2) * log(2*pi)
    res -= 0.5 * log_det_V
    res -= 0.5 * LinearAlgebra.dot(r, r)
    return(res)
end