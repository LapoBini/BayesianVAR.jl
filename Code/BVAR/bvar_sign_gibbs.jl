function bvar_sign_gibbs(
    Θ₁::Hyperparameter,           # Hyperpriors 
    data::DataFrame,              # data with dates and variable names
    yₛ::Array{Float64,2},         # dummy observations for dependent variables 
    xₛ::Array{Float64,2},         # dummy observations for independent variables
    Hᵢ::Int64,                    # horizon of the IRFs
    sign_s::Array{Float64,2},     # sign of the impulse 
    pos_policy::Vector{Float64};  # position of the shocked variables for each shock
    predictive_density = false,   # calculate predictive density 
    check_stationarity = true,    # test stationarity of the draw 
    max_draw           = 120000   # max. n. of draws of standart normal for sign id. 
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
    n_shock = size(sign_s,2);

    # Create Variables, intercept is last column of x
    x   = X_companion(array, Int(lags), intercept_c = true)' |> Array{Float64,2};
    y   = array[lags+1:end,:];
    t,k = size(y);
    tₐ  = size(yₛ,1);
    ν   = t + tₐ;

    # Allocate results
    PD  = zeros(H, k, reps-burnin);
    IRF = zeros(Hᵢ, k, reps-burnin, k);

    # Variables for the posteriors 
    Y = [y; yₛ];
    X = [x; xₛ]; 


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
    valid_draw = zeros(0);

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
                yₕ             = bvar_predictive_dens(k, H+1, t, lags, y, b, 𝛹)
                PD[:,:,jgibbs] = view(yₕ, p+1:H+p, :)
            end

            # ----------------------------------------------------------------------
            # 4 - Compute Matrix Contemporaneous Response
            # ----------------------------------------------------------------------
            rep_sign = 1;            # keep track number of repetitions of the draw
            ϕ        = zeros(k, k);  # pre-allocation of the matrix of t=0 responses
            rstr_met = 0;            # n. of shocks with all restr. met at each draw 
            flag     = 0;            # switching to 1 when all conditions are met

            # Cholesky decomposition of variance/covariance matrix 
            B₀ = Matrix(cholesky(Hermitian(𝛹)).U);

            # Loop until all restrictions are satisfied
            while flag < 1 && rep_sign <= max_draw

                # Draw a candidate matrix of contemporaneous restrictions 
                Q, R        = qr(randn(k,k)); # QR Factorization of a draw N(0, Iₖ)
                B̃₀          = Q * B₀;         # Matrix of contemporaneous shock Aˢ₀
                ϕ[1:k,1:k]  = B̃₀;            # Allocate the candidate  

                # Check if candidate matrix Aˢ₀ satisfies all the restrictions 
                for kk in 1:n_shock
 
                    # Position restriction and their signs
                    pos      = findall(pos_policy .== kk)[1];
                    rstr_pos = findall(.!isnan.(sign_s[:,kk]));
                    rstr_sgn = sign_s[rstr_pos,kk];

                    # Pick the restricted coefficient and n. of restrictions 
                    b₀ = B̃₀[pos, rstr_pos]
                    r  = length(rstr_pos)

                    # Check if the n. of satisfied restrictions equal the n. of
                    # restrictions specified for a given shock. The elseif check 
                    # if they meet the restrictions with the reversed signs.
                    if sum(sign.(b₀ .* rstr_sgn)) .== r
                        rstr_met += 1
                    elseif sum(sign.(b₀ .* rstr_sgn)) .== -r
                        rstr_met += 1
                        ϕ[1:k,1:k] = -B̃₀
                    end
                end

                # Check is all the shocks are satisfied, id est all the restrictions
                # of each shock are satisfied, otherwise keep iterating
                rstr_met == n_shock ? flag = 1 : rstr_met = 0; 

                # Update iterator for sign computation
                rep_sign += 1;

            end

            # ----------------------------------------------------------------------
            # 5 - Compute Structural IRFs by Recursion
            # ----------------------------------------------------------------------
            # Compute response for each column of the matrix of contemporaneous 
            # response B₀ (we will need for the Forecast Error Variance Decomposition)
            # but ONLY IF the draw of B₀ was succesful
            if flag == 1 
                for kk in 1:k

                    # select contemporaneous response to shock variable k 
                    B₀ = (ϕ[kk,:])' |> Array{Float64,2};

                    # Pre-allocate memory and create one unit shock, Hᵢ horizon IRFs
                    ŷ = zeros(lags+Hᵢ, k); # Save dynamic of the variables
                    ε = zeros(lags+Hᵢ, 1); # sequence of shocks (all zeros except h = 0)
                    ε[lags+1] = 1;         # shock to the policy measure at h = 0

                    # System is stationary, we consider deviation from steady state and
                    # before h = 0 the system was in steady state (variables were 0). For
                    # this reason we can also eliminate the intercept 
                    for tt = lags+1:Hᵢ+lags
                        x̂ = Array{Float64,2}(undef, lags, k)
                        for ji = 1:lags
                            x̂[ji,:] = ŷ[tt-ji,:]
                        end
                        x̂       = [reshape((x̂)', 1, k*lags) 0]
                        ŷ[tt,:] = x̂ * A + ε[tt] .* B₀;
                    end

                    # Allocate all the results
                    IRF[:,:,jgibbs,kk] = view(ŷ, lags+1:Hᵢ+lags, :)
                end
            
                # Valid draw index - repetition of post-burnin with valid B₀
                append!(valid_draw, jgibbs);

            end
        end
    end

    return PD, IRF, valid_draw

end
