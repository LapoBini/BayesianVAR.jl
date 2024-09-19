# ------------------------------------------------------------------------------
# HELPER FUNCTIONS FOR GIBBS SAMPLING AND STRUCTURAL IDENTIFICATION 
# ------------------------------------------------------------------------------
# These set of functions are used to implement gibbs sampling and identification
# via sign restrictions 
# Author: Lapo Bini lbini@ucsd.edu
# ------------------------------------------------------------------------------

function bvar_MNRD(
    βₐ::Array{Float64,1},   # posterior mean (constant)
    𝛹::Array{Float64,2},    
    Ωₓ⁻¹::Array{Float64,2},
    max_try::Int64,
    k::Int64,
    p::Int64;
    check_stationarity = true
    )

    # --------------------------------------------------------------------------
    # RANDOM DRAW MULTIVARIATE NORMAL AND CHECK STATIONARITY 
    # --------------------------------------------------------------------------

    # Pre-allocate space while loop
    problem   = 0;
    check     = -1;
    num_try   = 1;
    stability = 0;
    β⁺        = 0;

    # Compute variance of posterior distribution 
    V⁺ = kron(𝛹, Ωₓ⁻¹);

    while check < 0 && num_try < max_try

        # Random Draw from posterior distribution centered at posterior mean βₐ
        β⁺ = βₐ + (randn(1,k*(k*p+1))*Matrix(cholesky(Hermitian(V⁺)).U))';

        # Transform from vector to matrix to compute eigenvalues 
        kp    = k*p;
        kp₋₁  = k*(p-1);
        C_tmp = zeros(kp, kp);
        C_tmp[k+1:kp,1:kp₋₁] = eye(kp₋₁);

        # reshape
        tmp = reshape(β⁺, kp+1, k);
        tmp = view(tmp, 1:kp, 1:k)';
        C_tmp[1:k,1:kp] = tmp;

        # Eigenvalue matrix of coefficients
        λₛ = maximum(abs.(eigvals(C_tmp)));
        stability = λₛ >= 1

        # (i) Stationary BVAR: we are fine with the draw
        if stability == 0 || check_stationarity == false
            check = 10

        # (ii) Non-Stationary BVAR - keep iterating
        else
            num_try += 1
        end

    end

    # If after all the repetitions the system is non stable, report the issue
    # putting problem = 1
    if stability > 0 || check_stationarity == false
        problem = 1
    end

    return β⁺, problem

end


function bvar_predictive_dens(
    k::Int64,
    H::Int64,
    t::Int64,
    p::Int64,
    y::Array{Float64,2},
    b::Array{Float64,2},
    𝛹::Array{Float64,2}
    )

    # --------------------------------------------------------------------------
    # COMPUTE PREDICTIVE DENSITY 
    # --------------------------------------------------------------------------

    # Draw N(0,1) innovations for variance and mean equation:
    Σ = Matrix(cholesky(Hermitian(𝛹)).U)
    u = randn(H+p, k)

    # Note we only need H*K innovations, but adding an extra p draws makes
    # the indexing in the loop below a bit cleaner.
    yₕ        = zeros(H+p, k)
    yₕ[1:p,:] = y[t-p+1:t,:]

    # Generate forecasts
    for tt = p+1:H+p

        # Create Regressors
        xₕ = Array{Float64,2}(undef, p, k)
        for ji = 1:p
            xₕ[ji,:] = yₕ[tt-ji,:] 
        end

        # Prediction with innovations
        xₕ = [reshape((xₕ)', 1, k*p) 1]
        yₕ[tt,:] = xₕ * b + u[tt,:]' * Σ

    end

    return yₕ

end


function bvar_IWRD(
    S::Array{Float64,2},
    ν::Int64,
    )

    # --------------------------------------------------------------------------
    # RANDOM DRAW FROM INVERSE WISHART 
    # --------------------------------------------------------------------------
    # Draw from Inverted Wishert Distribution 𝛹 ∼ iW(S,ν) S is the (kxk) scale 
    # matrix and ν are the degrees of freedom. 𝛹 = (∑ₜ wᵢwᵢ')⁻¹ with wᵢ∼N(0,S)
    #
    # The Inverse Wishart distribution could be seen as multivariate version of 
    # the inverse gamma distribution (used for variance in univariate case) and
    #  it is a conjugate prior for the var/cov matrix of a multivariate normal.
    #
    # Random number generator for the inverse Wishart distribution): In order to
    # obtain a (nxn) random draw from an iW(S,ν)
    # 1 - Compute L the lower triangular Choleski of the scale matrix S = LL'.
    # 2 - Create Z = [z₁ … zᵥ] made by T draws of (kx1) random vectors zᵢ ∼ N(0,Iₖ)
    # 3 - Arrange vectors into a (Txn) matrix W = (L ⋅ Z)', where row i of W is  
    # 4 - estimate X = (W'⋅W)⁻¹ = (∑ₜ wᵢwᵢ')⁻¹ which is a random draw from iW(S,ν)
    # --------------------------------------------------------------------------

    # Allocation 
    k = size(S,1)
    W = zeros(ν,k)

    # 1 - Compute lower triangular matrix of Scale positive definite matrix ν
    L = Matrix(cholesky(Hermitian(S)).U)';

    # 2 & 3 - Draw independent normal vector zᵢ ∼ N(0,Iₖ) and create W = (L ⋅ Z)'
    for i = 1:ν
        W[i,:] = (L * randn(k,1))'
    end

    # 4 - Since (W'⋅W) ~ W(S,ν), now we want to obtain the inverse wishart:
    𝛹 = inv(W'*W);

    return 𝛹

end


function bvar_allign_shock(
    data::DataFrame, 
    instrument::DataFrame,
    lags::Int64
    )

    # --------------------------------------------------------------------------
    # ALLIGN DATES OF EXTERNAL INSTRUMENT AND VAR RESIDUALS
    # --------------------------------------------------------------------------
    # Dates of shock and var 
    date_z = instrument[:,1] .|> Date;
    ref_t  = data[lags+1:end,1] .|> Date;

    # Starting late, endind early
    beg         = [date_z[1]; ref_t[1]];
    fin         = [date_z[end]; ref_t[end]];
    start_shock = argmax(beg);
    end_shock   = argmin(fin);

    idx_u = findall((ref_t .>= beg[start_shock]) .& (ref_t .<= fin[end_shock]));
    idx_z = findall((date_z .>= beg[start_shock]) .& (date_z .<= fin[end_shock]));

    return idx_u, idx_z 

end


function bvar_irf_iv(
    pos_shock::Int64,
    ε::Matrix{Float64},
    b::Matrix{Float64},
    Z::Vector{Any},
    k::Int64,
    Tz::Int64,
    lags::Int64,
    Hᵢ::Int64
    )

    # --------------------------------------------------------------------------
    # STRUCTURAL IDENTIFICATION VIA EXTERNAL INSTRUMENT
    # --------------------------------------------------------------------------
    # Compute impulse response functions {Θ₀, ⋯ , Θₕ}, starting point is the 
    # Structural VAR and then we derive the reduce form:
    # B₀yₜ =  B₁yₜ₋₁ + … + Bₚyₜ₋ₚ + uₜ
    # yₜ   =  B₀⁻¹B₁yₜ₋₁ + … + B₀⁻¹Bₚyₜ₋ₚ + B₀⁻¹uₜ
    # yₜ   =  Φ₁yₜ₋₁ + … + Φₚyₜ₋ₚ + εₜ
    # Relation between reduced form error and structural shocks: εₜ = B₀⁻¹uₜ
    # Identification procedure: by instrumental variable, ory Proxy-SVAR. Wold
    # representation of the Companion form VAR:
    # Yₜ   =  Ω₁Yₜ₋₁ + Eₜ
    # Yₜ   = ∑ Ωʰ Eₜ
    # J Yₜ = ∑ J Ωʰ J'J Eₜ
    # yₜ   = ∑ 𝛹ₕ εₜ =  ∑ 𝛹ₕ B₀⁻¹ B₀ εₜ
    # yₜ   = ∑ Θₕ uₜ
    # --------------------------------------------------------------------------
    # Construct Dynamic Multiplier and pre allocate B₀⁻¹
    Φ    = [b[1:end-1,:]'; eye(k*(lags-1)) zeros(k*(lags-1), k)];
    B₀⁻¹ = zeros(k, k);

    # Estimate reduced form residual and first stage estimation B₀⁻¹
    E_εZ  = ε' * Z;
    E_ε₁Z = E_εZ[pos_shock:pos_shock,:];      # relevance condition
    bⱼ    = (E_εZ * E_ε₁Z)./(E_ε₁Z * E_ε₁Z);  # first stage coefficient

    # Allocate contemporaneous response in the column of the instrumented variable
    B₀⁻¹[:,pos_shock] = bⱼ;

    # Selection matrix and allocate memory for results
    J   = [eye(k) zeros(k, k*(lags-1))]; # selection matrix
    IRF = zeros(Hᵢ, k);

    # Estimation Impulse response function
    IRF[1,:] = B₀⁻¹[:,pos_shock]' |> any2float;
    for h in 1:Hᵢ-1
        Ψ   = J * Φ^h * J'
        irf = (Ψ * B₀⁻¹)[:,pos_shock]
        IRF[h+1,:] = irf;
    end

    return IRF

end