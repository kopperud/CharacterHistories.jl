export Brownian

struct Brownian
    sigma2::Float64
    mean::Float64
    observation_variance::Float64
end

Brownian(sigma2, mean) = Brownian(sigma2, mean, 0.0)