function bvar_marginal_likelihood(
    y::Array,  # oberved LHS
    x::Array,  # observed RHS
    yₛ::Array, # dummy observation LHS
    xₛ::Array  # dummy observation RHS
    )

    # ----------------------------------------------------------------------
    # COMPUTE MARGINAL LOG LIKELIHOOD
    # ----------------------------------------------------------------------
    # Let's called the model parameters ϕ = {β, Ψ}, from the Bayes'Theorem:
    #
    #                  π(ϕ|Y) = [f(Y|ϕ) ⋅ π(ϕ)] ⋅ m(Y)⁻¹ 
    #          log(π(ϕ|Y)) = log(f(Y|ϕ)) + log(π(ϕ)) - log(m(Y))
    #          log(m(Y)) = log(f(Y|ϕ)) + log(π(ϕ)) - log(π(ϕ|Y))
    #   log(m(Y)) = log(f(Y|ϕ)) + log(π(ϕ)) - log(π(β|Ψ,Y)) - log(π(Ψ|Y))
    #
    # but m(Y) = m(yₛ,y) since Y is made by dummy observation (yₛ) and the 
    # observed data (y). Moreover, yₛ and y are independent, ie 
    # m(Y) = m(yₛ,y) = m(yₛ)⋅m(y). Therefore:
    #
    #                           m(y) = m(Y)/m(yₛ)
    #
    # given the sum of coefficient and the dummy initial observation priors
    # the marginal likelihood should be computed on the original data only. 
    # As derived in Giannone, Lenza, and Primiceri (2014), this is equivalent 
    # to taking the ratio between the marginal likelihood of the augmented 
    # data set relative to the marginal likelihood of the dummy observations
    # ----------------------------------------------------------------------

    # Sample size (time t and variables k) and dof prior and posterior
    t,k = size(y)
    v₀  = size(yₛ,1) 
    v₁  = v₀ + t    

    # Allocate memories
    β₀ = [];
    Ω₀ = [];

    # Prior Moments (scale matrix Σ₀ for iW(Σ₀, v₀) of prior distributions)
    try β₀ = inv((xₛ)'*xₛ)*((xₛ)'*yₛ) catch e β₀ = pinv((xₛ)'*xₛ)*((xₛ)'*yₛ) end
    ϵ₀ = yₛ - xₛ * β₀
    Σ₀ = ϵ₀'*ϵ₀ 

    # Variance Covariance Matrix
    try Ω₀ = inv((xₛ)'*xₛ) catch e Ω₀ = pinv((xₛ)'*xₛ) end
    P = inv(eye(t) + x * Ω₀ * x')

    # Compute the log of the ratio mult. Gamma function of prior and posterior 
    lnΓ = log_multivariate_gamma(k,v₀) - log_multivariate_gamma(k,v₁)

    # compute log(m(y)) = log(m(Y)) - log(m(yₛ))
    Π   = log_marginal_likelihood(lnΓ,t,k,v₀,v₁,P,Σ₀,y,x,β₀)

    return Π

end
