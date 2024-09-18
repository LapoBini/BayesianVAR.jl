function bvar_opt_prior(
    Θ::Hyperparameter,
    Xm::Array{Float64,2},
    μ::Array{Float64,1},
    σ::Array{Float64,2},
    δ::Array{Float64,2}
    )

    # --------------------------------------------------------------------------
    # OPTIMIZATION PRIOR TIGHTNESS 
    # --------------------------------------------------------------------------
    # The functions is used to compute the optimal informativeness of the prior
    # based on the marginal log-likelihood. The prior are parametrized using
    # the hyperpriors parameters which are treated as additional parameters, in 
    # the spirit of hierarchical modeling. No estimated posteriors are used to 
    # compute the marginal likelihood.
    #
    # Priors: vec(β)|Ψ ∼ N(vec(β₀), Ψ ⊗ Ω₀) and Ψ ∼ iW(Σ₀,v₀) are 
    # implemented using dummy observations, where β₀ is computed by OLS on 
    # dummy observation, Σ₀ from the corresponding residuals and Ω₀ = (xₛ'xₛ)⁻¹
    #
    # Posteriors: vec(β)|Ψ ∼ N(vec(β̃), Ψ ⊗ Ω̃) and Ψ ∼ iW(Σ₁,v₁) where
    # where β̃ is computed by OLS on the dummy augmented sample, Σ₁ from the 
    # corresponding residuals and Ω̃ = (x̃'x̃)⁻¹
    #
    # Author: Lapo Bini, lbini@ucsd.edu
    # --------------------------------------------------------------------------

    @unpack λ, τ, γ, ε, p, H, reps, burnin, max_try, update = Θ;
    t,k = size(Xm);
    ℒ   = zeros(length(p),length(λ));

    # Allocate arrays
    opt_p  = zeros(length(p),length(λ));
    opt_λ  = copy(opt_p);
    opt_τ  = copy(opt_p);
    opt_γ  = copy(opt_p);

    for i = 1:length(p)

        # Create the arrays for the estimation  of the likelihhod
        x  = X_companion(Xm, p[i], intercept_c = true)' |> Array{Float64,2};
        y  = Xm[p[i]+1:end,:];

        # Compute value marginal likelihood for different value prior
        for j = 1:length(λ)

            # Create dummies and compute the marginal likelihood
            yₛ, xₛ = bvar_dummies(λ[j], τ[j], γ[j], ε, p[i], k, δ, μ, σ);
            Π     = bvar_marginal_likelihood(y, x, yₛ, xₛ);

            # Allocate results
            ℒ[i,j]     = Π
            opt_p[i,j] = p[i]
            opt_λ[i,j] = λ[j]
            opt_τ[i,j] = τ[j]
            opt_γ[i,j] = γ[j]
        end

    end

    # Select the best hyperprior specification based on the marginal likelihood
    # We do a grid-search 
    id = findmax(ℒ)[2]
    p̂  = opt_p[id]
    λ̂  = opt_λ[id]
    τ̂  = opt_τ[id]
    γ̂  = opt_γ[id]

    return p̂, λ̂, τ̂, γ̂

end
