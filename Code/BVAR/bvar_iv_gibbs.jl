function bvar_iv_gibbs(
    Θ₁::Hyperparameter,           # Hyperpriors 
    data::DataFrame,              # data with dates and variable names
    yₛ::Array{Float64,2},         # dummy observations for dependent variables 
    xₛ::Array{Float64,2},         # dummy observations for independent variables
    Hᵢ::Int64,                    # horizon of the IRFs
    instrument::DataFrame,        # date of the instrument and the shock
    tickers::Vector{Any};         # ticker to identify the column for the shock
    predictive_density = false,   # calculate predictive density 
    check_stationarity = true,    # test stationarity of the draw 
    )

    # ------------------------------------------------------------------------------
    # COMPUTE IRF SIGN RESTRICTIONS, Author: Lapo Bini, lbini@ucsd.edu
    # ------------------------------------------------------------------------------

    # Unpack and auxiliary variables 
    @unpack λ, τ, γ, ε, p, H, reps, burnin, max_try, update = Θ₁
    lags    = p[1];
    array   = data[:,2:end] |> any2float;
    names_v = names(data)[2:end];
    dates   = data[lags+1:end,1] .|> Date;

    # Create Variables, intercept is last column of x
    x   = X_companion(array, Int(lags), intercept_c = true)' |> Array{Float64,2};
    y   = array[lags+1:end,:];
    t,k = size(y);
    tₐ  = size(yₛ,1);
    ν   = t + tₐ;

    # Allocate results
    PD  = zeros(H, k, reps-burnin);
    IRF = zeros(Hᵢ, k, reps-burnin, 1);

    # Variables for the posteriors 
    Y = [y; yₛ];
    X = [x; xₛ]; 


    # ------------------------------------------------------------------------------
    # ALLIGN SHOCK AND DATA
    # ------------------------------------------------------------------------------
    idx_u, idx_z = bvar_allign_shock(data, instrument, lags);
    Z            = instrument[idx_z,2];
    pos_shock    = findall(tickers .== names(instrument)[2])[1];
    Tz           = length(Z);


    # ------------------------------------------------------------------------------
    # POSTERIOR DISTRIBUTION
    # ------------------------------------------------------------------------------
    Ωₓ⁻¹ = inv(X'*X);
    A    = Ωₓ⁻¹*(X'*Y); # Reduced form coefficient (mean of posterior distribution)
    u    = Y - X*A;     # residual using posterior mean of coefficients 
    𝛹    = (u' * u)./ν; # Scale matrix for Inverse wishart (it will be updated)

    # Iterators, progress and index of iterations post-burnin with valid B₀ draws
    jgibbs     = 0;
    repetition = [1; collect(1:1:100) .* (reps/100) |> Array{Int64,1}];
    percentage = [0; collect(1:1:100)];

    # Auxiliary variables for jibb sampling 
    βₐ = deepcopy(A[:]); # fixed parameter (posterior mean of coefficient)
    β₁ = deepcopy(A[:]); # changing one (it will be updated with random draws)
    B₀ = [];             # used for sign restriction 


    # ------------------------------------------------------------------------------
    # GIBBS SAMPLING
    # ------------------------------------------------------------------------------
    for i in 1:reps

        # Display progress:
        idx = findall(repetition .== i);
        if ~ isempty(idx)
            ite = percentage[idx[1]]
            i <= burnin ? ite_aux = "Burnin" : ite_aux = "Post-burnin";
            println("BVAR > Structural Analysis > Gibbs Sampling: $ite_aux > $ite% iteration")
        end

        # --------------------------------------------------------------------------
        # 1 - Random Draw BVAR Coefficients and Check Stationarity 
        # --------------------------------------------------------------------------
        β, problem = bvar_MNRD(βₐ, 𝛹, Ωₓ⁻¹, max_try, k, lags);
        problem == 1 ? β = βₐ : β₁ = β; # problem = 0 BVAR stationary, draw accepted

        # --------------------------------------------------------------------------
        # 2 - Random Draw Variance Covariance Matrix
        # --------------------------------------------------------------------------
        # We take the inverse of Var/Cov of the residual because we draw initially 
        # from the wishart distribution with scale parameter the inverse of the 
        # original scale matrix for IW, to then take the inverse of the draw. 
        b  = reshape(β₁, k*lags+1, k); # Coefficient after draw matrix form
        uⱼ = Y - X*b;                 # residuals 
        S  = inv(uⱼ' * uⱼ);           # Inverse of Var/Cov matrix residual
        𝛹  = bvar_IWRD(S, ν);         # S = scale matrix, ν = degrees of freedom 

        # Sometimes the scale matrix is defined as S = Σ₀ + (Yₜ - Xₜ⋅A₁)'(Yₜ - Xₜ⋅A₁)
        # but here we are using dummy variables, so all together and not the sum of 
        # prior variance plus var/cov of reduced form residuals:
        # S = (Y⁺ - X⁺A₁)'(Y⁺ - X⁺A₁) where Y⁺ and X⁺ appended data


        # --------------------------------------------------------------------------
        # POST BURNIN PERIOD 
        # --------------------------------------------------------------------------
        if i > burnin && problem == 0
            jgibbs += 1

            # ----------------------------------------------------------------------
            # 3 - Compute predictive density 
            # ----------------------------------------------------------------------
            if predictive_density
                yₕ             = bvar_predictive_dens(k, H+1, t, lag, y, b, 𝛹)
                PD[:,:,jgibbs] = view(yₕ, p+1:H+p, :)
            end

            # ----------------------------------------------------------------------
            # 4 - Compute Structural IRF via External Instrument
            # ----------------------------------------------------------------------
            # Take residual corresponding to the same period of the shock
            ε = uⱼ[idx_u,:];

            # Compute IRF 
            IRF[:,:,jgibbs, 1] = bvar_irf_iv(pos_shock, ε, b, Z, k, Tz, lags, Hᵢ);
        end
    end

    return PD, IRF, pos_shock

end
