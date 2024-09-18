# ------------------------------------------------------------------------------
# HELPERS FOR GAMMA FUNCTION
# ------------------------------------------------------------------------------
# These set of functions are used to optimization of hyperprior and to draw from
# the Inverse Wishart Distribution which is a conjugate prior from the normal 
# distribution.
# Author: Lapo Bini lbini@ucsd.edu
# ------------------------------------------------------------------------------

function sum_log_gamma(
    k::Int64, # number of random variable 
    v::Int64  # degrees of freedom 
    )

    # ------------------------------------------------------------------------------
    # Compute sum of log(Γₖ(v)) with Γₖ(v) multivariate gamma function with K equals 
    # to the number of variables in the model.
    #
    # Give, log(Γₖ(v)) = k(k-1)⋅log(π) + ∑ᵢᵏ log(Γ(v+(1-i)/2)) where Γ(v+(1-i)/2) is
    # the univariate gamma function. This function compute the sum over the 
    # univariate log-gamma γ = ∑ᵢᵏ log(Γ(v+(1-i)/2))
    # ------------------------------------------------------------------------------

    # Pre-allocate value 
    γ = 0;

    # Compute value of the sum of the log univariate gamma 
    for i = 1:k
        γ += logabsgamma((v+1-i)/2)[1]
    end

    return γ

end


function log_multivariate_gamma(
    k::Int64,
    v::Int64
    )

    # --------------------------------------------------------------------------
    # Compute the log of the multivariate (of order k) gamma function: log Γₖ(v)
    # --------------------------------------------------------------------------

    # log(Γₖ(v)) = c + γ where c = k(k-1)⋅log(π) and γ = ∑ᵢᵏ log(Γ(v+(1-i)/2))
    c  = (k*(k-1)/4)*log(pi)
    γ  = sum_log_gamma(k, v)
    c += γ

    return c

end


function log_marginal_likelihood(
    lnΓ::Float64,         # ratio of log multivariate gamma function
    t::Int64,             # number of observation 
    k::Int64,             # dimension of the model (n. of variables)
    v₀::Int64,            # dof posterior distribution iW(Σ₀,v₀)
    v₁::Int64,            # dof posterior distribution iW(Σ₁,v₁)
    P::Array{Float64,2},  # P = (Iₖ + xΩ₀x') with |P| = |Ω₀||Ω̃|⁻¹
    Σ₀::Array{Float64,2}, # Scale parameter for prior iW(Σ₀,v₀)
    y::Array{Float64,2},  # observed data (no dummy) - LHS
    x::Array{Float64,2},  # observed data (no dummy) - RHS
    β₀::Array{Float64,2}, # prior mean 
    )

    # --------------------------------------------------------------------------
    # Compute the marginal log likelihood following Giannone et Al. (2014)
    # For the derivation look at BEAR (2016) pag.241 Appendix A6. 
    # The marginal likelihood of the observed data is:
    #
    # log(m(y)) = log(m(Y))/log(m(yₛ)) 
    # log(m(y)) = - K + 0.5 log(|Ω₀||Ω̃|⁻¹) + 0.5 (v₀/v₁) log(|Σ₀||Σ₁|⁻¹)
    # --------------------------------------------------------------------------
    
    # Normalizing constant: K = log(Γₙ(v₁/2)) - log(Γₙ(v₀/2)) - (tK/2)⋅log(π)
    K  = (lnΓ+(t*k/2)*log(pi))

    # Scale parameter posterior distribution 
    Σ₁ = Σ₀+(y-x*β₀)'*P*(y-x*β₀)
    
    # Difference marginal loglikelihood stacked and dummy observations 
    Π  = -K + 0.5*k*log(det(P)) + (v₀/2)*log(det(Σ₀)) - (v₁/2)*log(det(Σ₁))

    return Π
end

