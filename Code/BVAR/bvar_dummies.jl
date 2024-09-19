function bvar_dummies(
    λ::Float64,
    τ::Float64,
    γ::Float64,
    ε::Float64,
    p::Int64,
    k::Int64,
    δ::Array{Float64,2},
    μ::Array{Float64,1},
    σ::Array{Float64,2}
    )

    # ------------------------------------------------------------------------------
    # CREATE DUMMY OBSERVATIONS TO IMPLEMENT PRIORS
    # ------------------------------------------------------------------------------
    # Prior distributions are vec(β)|Ψ ∼ N(vec(β₀), Ψ ⊗ Ω₀⁻¹) and Ψ ∼ iW(Σ₀,v₀). This 
    # function create the dummy variables to construct the prior distributions. It 
    # makes use of the parameters related to the Minnesota prior estimated before. 
    #
    # 1 - Minnesota Prior: β₀ = [A₁ A₂ … Aₚ] are assumed to be a priori 
    #     independent and normally distributed where E[(Aₖ)Iᵢⱼ] = δᵢ if j=i and k=1 
    #     and E[(Aₖ)Iᵢⱼ] = 0 otherwise. δᵢ = 1 or 0 depending on mean reversion. 
    #
    # 2 - Prior on the Intercept c is diffuse (governed by a veru small number ε)
    #
    # 3 - Prior on the variance: Ψ = diag(σ₁²,…,σₖ²) and Ω₀ = (Xₛ'Xₛ) is chosen so 
    #     that the  Ψ ⊗ Ω₀⁻¹ is made by Var((Aₖ)Iᵢⱼ) = λ²/k² if i=j and
    #     Var((Aₖ)Iᵢⱼ) = ϑ⋅(λ²/k²)⋅(σᵢ²/σⱼ²) otherwise
    #
    # 4 - Sum-of-Coefficients prior: it implements the inexact differencing, id est
    #     ΔYₜ = c - (Iₙ - A₁ - A₂ - … - Aₚ) Y₍₋₁ + B₁ ΔYₜ₋ₗ + … + Bₚ ΔYₜ₋ₚ + Δuₜ
    #     with (Iₙ - A₁ - A₂ - … - Aₚ) = 0 is implied by a VAR in first difference. 
    #     By inxact differncing we are going to shrink this value toward zero. 
    #     A shortcoming of the sum-of-coefficients strategy is that it rules out 
    #     cointegration in the limit, which may be undesirable:
    #
    #      ΔYₜ = c - (Iₙ - A₁ - A₂ - … - Aₚ) Yₜ₋₁ + B₁ ΔYₜ₋ₗ + … + Bₚ ΔYₜ₋ₚ + Δuₜ
    #                    ΔYₜ = c + B₁ ΔYₜ₋ₗ + … + Bₚ ΔYₜ₋ₚ + Δuₜ
    #                 Yₜ = Yₜ₋₁ + c + B₁ ΔYₜ₋ₗ + … + Bₚ ΔYₜ₋ₚ + Δuₜ
    #
    #     All yₜ have unit root but withouth an error correction term we cannot have 
    #     cointegration, that is why we will use the next prior. 
    #
    # 5 - Dummy initial observation prior: shrinks the forecast of each variable at
    #     the beginning of the sample toward a no-change forecast. This prior 
    #     solves the cointegration issue of the sum-of-coefficient prior. 
    #
    # Parameter λ controls the overall tightness of the prior distribution around
    # the random walk or white noise and governs the relative importance of the 
    # prior beliefs with respect to the information contained in the data. λ = 0
    # the posterior is equal to the prior, while λ = ∞ the posterior is the OLS.
    # This parameter should increase with the number of variables included.
    # 
    # The factor 1/k² is the rate at which prior variance decreases with increasing 
    # lag length of the VAR
    #
    # (σᵢ²/σⱼ²) accounts for the different scale and variability of the data.
    #
    # ϑ ∈ [0,1] governs the extent to which the lags of other variables are “less 
    # important” than the own lags. In the context of the SVAR  we need to take into 
    # account possible correlation among the residual of different variables. 
    # Consequently, Litterman’s assumption of fixed and diagonal covariance matrix
    # is problematic. To overcome this problem we follow Kadiyala and Karlsson (1997) 
    # and Robertson and Tallman (1999) and impose a Normal inverted Wishart prior 
    # which retains the principles of the Minnesota prior. This is possible under 
    # the condition that ϑ = 1 which is assumed in what follows.
    #
    # τ = 10⋅λ controls for the degree of shrinkage, as τ goes to zero we approach
    # the case of exact differences and, as τ goes to ∞ , we approach the case of no
    # shrinkage.
    #
    # γ controls the tightness of the dummy initial observation prior, when γ 
    # approach infinity we have no shrinkage
    # ------------------------------------------------------------------------------

    Jₚ = diagm(1:p);
    
    # 1 - Minnesota Prior 
    # Auxiliary variable LHS
    yₐ = [(σ.*δ)./λ; zeros(k*(p-1),k);   # Restriction on autoregressive coefficients
          σ;                             # prior for the covariance matrix
          zeros(1,k)];                   # Uninformative prior for intercept
    
    # Auxiliary variable RHS
    xₐ = [kron(Jₚ, σ./λ) zeros((k*p),1); # Restriction on autoregressive coefficients
          zeros(k,(k*p)+1);              # prior for the covariance matrix
          zeros(1,k*p) ε];               # Uninformative prior for intercept

    # 2 - Sum-of-Coefficient Prior
    # (Equation in Banbura et al. is wrong,it's not (1 2 … p) but a vector 1̲ₚ that 
    # enters the kronecker product for xₐ₂)
    yₐ₂ = δ.*μ./τ;
    xₐ₂ = [kron(ones(1,p),yₐ₂) zeros(k,1)];

    # 3 - Dummy-Initial-Observation Prior 
    yₐ₃ = μ'./γ;
    xₐ₃ = [kron(ones(1,p),yₐ₃) 1 ./γ];

    # Stack together 
    yₐ = vcat(yₐ,yₐ₂,yₐ₃)
    xₐ = vcat(xₐ,xₐ₂,xₐ₃)

    return yₐ, xₐ

end
