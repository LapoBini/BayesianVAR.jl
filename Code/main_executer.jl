# ------------------------------------------------------------------------------
# EXECUTER BVAR
# ------------------------------------------------------------------------------
# Set-up of the model 
include(pwd()*"/Code/SetUp/SetupBVAR.jl");

# run the model 
bvar_main(data_path, start_date, end_date, results_folder,
# If you want set manually the hyperpriors 
          λ = [0.1,1.0,10.0,100],     # shrinkage minnesota prior
          τ = [0.1,1.0,10.0,100].*10, # shrinkage sum-of-coefficients 
          γ = [0.1,1.0,10.0,100].*10, # shrinkage dummy-initial-observations 
          ε = 0.0001,                 # diffuse prior intercept
          p = [2,3,4,5,6,7,8,12]      # lag-length BVAR 
          )