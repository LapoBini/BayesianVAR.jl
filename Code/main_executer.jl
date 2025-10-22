# ------------------------------------------------------------------------------
# EXECUTER BVAR
# ------------------------------------------------------------------------------
# Set-up of the model 
include(pwd()*"/Code/SetUp/SetupBVAR.jl");

# run the model 
bvar_main(data_path, start_date, end_date, results_folder,
          # If you want set manually the hyperpriors 
          λ = [0.1,1.0,10.0,100],     # shrinkage minnesota prior
          τ = [0.1,1.0,10.0,100].*10, # shrinkage sum-of-coefficients ptiot 
          γ = [0.1,1.0,10.0,100].*10, # shrinkage dummy-initial-observations prior 
          ε = 0.0001,                 # diffuse prior intercept
          p = [4],                    # lag-length BVAR 
          # length impulse response functions 
          years = 5,
          # If you want historical decomposition with sign restrictions
          historical_decomp = ["GDP"; "INFL"] # (put seriesID from Legend table)
)