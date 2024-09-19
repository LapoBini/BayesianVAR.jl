function bvar_main(
    data_path::String,
    start_date::String,
    end_date::String,
    results_folder;
    # Hyperpriors for, H is horizon predictive density 
    λ       = [0.1,1.0,10.0,100],
    τ       = λ.*10,
    γ       = λ.*10,
    ε       = 0.0001,
    p       = [2,3,4,5,6,7,8,12],
    reps    = 10000,
    burnin  = 3000,
    max_try = 2000,
    update  = 1000,
    # Extra
    predictive_density = false,
    check_stationarity = true,
    # Event study with external instrument shock:
    name_iv_shock   = "Oil",
    forecast_origin = []     # ["30/06/2021"; "31/01/2021"]
    )

    # --------------------------------------------------------------------------
    # BAYESIAN VAR EXECUTER 
    # --------------------------------------------------------------------------
    # The function implements the BVAR which consists if a combination of the 
    # most commonly used priors in the literature: MINNESOTA, SUM-OF-COEFFICIENTS 
    # (Doean et al., 1984) and DUMMY INITIAL OBSERVATIONS (Sims, 1992) priors. 
    # The priors are implemented through dummy observations. The priors are 
    # consistent with unit root and cointegration processes.
    #
    # The prior distributions are implemented using Normal - Inverted Wishart
    # distributions, which are conjugate priors:
    #
    #             vec(β)|Ψ ∼ N(vec(β₀), Ψ ⊗ Ω₀⁻¹) and Ψ ∼ iW(Σ₀,v₀) 
    #
    # The Bayesian VAR is related to the one described in Banbura et al. (2007),
    # with the only difference that they didn't implement dummy initial obs.
    # The choice of the informativeness of the prior is based on Giannone, Lenza,
    # Primiceri (2012): prior distribution πᵧ(ϕ) with ϕ = {β, Ψ} vector of model
    # parameters and γ vector of hyperparameters, i.e. those coefficients that
    # parameterize the prior distribution, but do not directly affect the 
    # likelihood.
    #
    # Author: Lapo Bini, lbini@ucsd.edu
    # --------------------------------------------------------------------------

    # Recession dates and sample dates
    start_s = DateTime(start_date, "dd/mm/yyyy");
    end_s   = DateTime(end_date, "dd/mm/yyyy");


    # --------------------------------------------------------------------------
    # 1 - Load Dataset 
    # --------------------------------------------------------------------------
    println("BVAR > Read Data > Loading Dataset")
    data, ref_dates, tickers, prior, sign_s, name_s, instrument, pos_policy, 
    base_frq = readdata_haver(data_path, start_date, end_date);

    # Horizon for Impulse Response Functions (+1 because shock at time 0)
    base_frq == "m" ? Hᵢ = 61 : Hᵢ = 21;
    H = Hᵢ - 1;


    # --------------------------------------------------------------------------
    # 2 - Optimization Hyperpriors 
    # --------------------------------------------------------------------------
    # We are going to chose the hyperpriors (or hyperparameters) which maximize 
    # the marginal likelihood of the data given the priors. 
    println("BVAR > Hyperpriors > Optimization marginal likelihood")
    Θ  = Hyperparameter(λ, τ, γ, ε, p, H, reps, burnin, max_try, update);
    Xm = data[:,2:end] |> any2float;

    # Find the optimal parametrization of prior distribution using marginal 
    # distribution of the data 
    t, k       = size(Xm);
    μ, σ, δ, ϵ = bvar_prior(Xm, prior, 1);
    p, λ, τ, γ = bvar_opt_prior(Θ, Xm, μ, σ, δ);
    println("BVAR > Hyperpriors > Results > p = $p / λ = $λ / τ = $τ / γ = $γ")

    # Create the dummy observation for the optimal priors
    println("BVAR > creation dummy observations")
    Θ₁    = Hyperparameter([λ], [τ], [γ], ε, [p], H, reps, burnin, max_try, update);
    yₛ, xₛ = bvar_dummies(λ, τ, γ, Θ.ε, Int64(p), k, δ, μ, σ);


    # --------------------------------------------------------------------------
    # 3 - Structural Analysis: Sign Restrictions
    # --------------------------------------------------------------------------
    # Structural analysis by sign restriction that must be specified in the 
    # excel file of the data. Create a worksheet "Sign", first two
    # columns must contain tickers, name variables, and from the third on 
    # all the different signs (first row must have name of the shock)
    if ~isempty(sign_s)
    
        # Create Results folder 
        ind_dir   = readdir("./");
        "Results" in ind_dir ? nothing : mkdir("./Results");
        ind_dir   = readdir("./Results");
        "$results_folder" in ind_dir ? nothing : mkdir("./Results/$results_folder");

        # Create Folder for IRFs and FEVDs sign restrictions 
        list_dir = readdir("./Results/$results_folder");
        res_path = "./Results/$results_folder/IRF_sign";
        if size(findall(list_dir.==["IRF_sign"]),1) == 0
            mkdir(res_path);
        end

        # Estimation
        println("BVAR > Structural Analysis > Sign Identification")
        PD, IRF, valid_draw = bvar_sign_gibbs(Θ₁, data, yₛ, xₛ, Hᵢ, sign_s, pos_policy,
                                              check_stationarity = check_stationarity,
                                              predictive_density = predictive_density);

        # Documentation
        println("BVAR > Structural Analysis > Plot results")
        bvar_documentation(IRF, data, H, name_s, pos_policy, base_frq, res_path);
    end


    # --------------------------------------------------------------------------
    # 4 - Structural Analysis: External Instrument
    # --------------------------------------------------------------------------
    # Estimate the BVAR, compute the residuals and then estimate the desired
    # column of B₀⁻¹ corresponding to the instrumented variable using moment
    # condition  
    if ~isempty(instrument)

        # Create Results folder 
        ind_dir   = readdir("./");
        "Results" in ind_dir ? nothing : mkdir("./Results");
        ind_dir   = readdir("./Results");
        "$results_folder" in ind_dir ? nothing : mkdir("./Results/$results_folder");

        # Create Folder for IRFs and FEVDs sign restrictions 
        list_dir = readdir("./Results/$results_folder");
        res_path = "./Results/$results_folder/IRF_iv";
        if size(findall(list_dir.==["IRF_iv"]),1) == 0
            mkdir(res_path);
        end

        # Compute Structural BVAR and IRF using external Instrument 
        println("BVAR > Structural Analysis > IV Identification")
        PD, IRF, pos_shock = bvar_iv_gibbs(Θ₁, data, yₛ, xₛ, Hᵢ, instrument, tickers,
                                           check_stationarity = check_stationarity,
                                           predictive_density = predictive_density);

        # Documentation
        println("BVAR > Structural Analysis > Plot results")
        bvar_documentation(IRF, data, H, [name_iv_shock], pos_shock, base_frq, res_path);
    end
end


# --------------------------------------------------------------------------
# MANUAL RUN VARIABLES:
#=--------------------------------------------------------------------------
λ=[0.1, 1, 10.0,100.0]; τ=λ.*10; γ=λ.*10; ε=0.0001; p=[2,3,4,5,6,7,8,12]; H=60; 
burnin=2000; max_try=100; update=1000; reis=false; 
reps=5000; save_transformed=false; sign_restriction=false; Hᵢ=12*5;
predictive_density=false; check_stationarity=true; frequency = "Month";
forecat_origin = []; name_iv_shock = "Oil";
# ------------------------------------------------------------------------=# 
