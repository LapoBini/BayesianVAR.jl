function bvar_prior(
    y::Array{Float64,2},     # dataset 
    prior::Array{Float64,1}, # prior: 1 random wal, 0 white noise 
    p::Int64                 # lag-length: used to compute unconditional mean 
    )                

    # ------------------------------------------------------------------------------
    # MINNESOTA PRIOR IMPLEMENTATION: SHRINKAGE TOWARDS RW OR WN
    # ------------------------------------------------------------------------------
    # Prior distributions are vec(β)|Ψ ∼ N(vec(β₀), Ψ ⊗ Ω₀) and Ψ ∼ iW(Σ₀,v₀), this 
    # function create all the relevant variables to implement the Minnesota priors,
    # id est shrinkage towards random walk or white noise depending on mean 
    # reversion. All the variables will be use afterwards to create dummy variables.
    #
    # 1 - Mean prior of coefficients:
    #     β₀ = [A₁ A₂ … Aₚ] are assumed to be a priori independent and normally 
    #     distributed where E[(Aₖ)Iᵢⱼ] = δᵢ if j=1 k=1 and E[(Aₖ)Iᵢⱼ] = 0 
    #     otherwise in order to implement random walk Yₜ = c + δ Yₜ₋₁ + ε or white 
    #     noise Yₜ = c + ε
    # 2 - Prior on the intercept c is diffuse.
    # 3 - Prior on the variance: Ψ = diag(σ₁²,…,σₖ²)
    #
    # Author: Lapo Bini, lbini@ucsd.edu
    # ------------------------------------------------------------------------------

    t,k = size(y);
    μ   = zeros(k);
    σ   = zeros(k, k);
    δ   = diagm(prior);
    ϵ   = Array{Float64,2}(undef, t-p, k);

    x = X_companion(y, 1, intercept_c = false);
    Ι = ones(t-p);

    for i = 1:k
        yₜ = y[p+1:end,i];
        xₜ = x[i,:];
        βₖ = δ[i,i];

        μ[i]   = mean(yₜ - βₖ .* xₜ);
        X      = [xₜ Ι];
        ε      = yₜ - X * [βₖ, μ[i]];
        σ[i,i] = ε' * ε /length(yₜ);
        ϵ[:,i] = ε;
    end

    return μ, σ, δ, ϵ

end
