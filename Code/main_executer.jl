# ------------------------------------------------------------------------------
# EXECUTER BVAR
# ------------------------------------------------------------------------------
include(pwd()*"/Code/SetUp/SetupBVAR.jl");
bvar_main(data_path, start_date, end_date, results_folder)