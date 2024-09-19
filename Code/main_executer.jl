# ------------------------------------------------------------------------------
# EXECUTER BVAR
# ------------------------------------------------------------------------------
# Set-up of the model 
include(pwd()*"/Code/SetUp/SetupBVAR.jl");

# run the model 
bvar_main(data_path, start_date, end_date, results_folder,
# If you want set manually the hyperpriors 
λ  = [0.1,1.0,10.0,100],
τ  = [0.1,1.0,10.0,100].*10,
γ  = [0.1,1.0,10.0,100].*10,
ε  = 0.0001,
p  = [2,3,4,5,6,7,8,12]
)