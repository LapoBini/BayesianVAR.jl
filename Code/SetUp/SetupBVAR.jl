# ------------------------------------------------------------------------------
# 1 - Load Packages
# ------------------------------------------------------------------------------
println("BVAR Set-up > Load packages")
using DataFrames, LinearAlgebra, Dates, Statistics, Plots, LaTeXStrings, GLM
using XLSX, CSV, JLD2, FredData, Statistics, Dierckx, SpecialFunctions, Random
using GLM, Parameters 
Random.seed!(1234)


# ------------------------------------------------------------------------------
# 2 - Setup BVAR
# ------------------------------------------------------------------------------
println("BVAR Set-up > Manual settings")

# Excel file to load
data_file = "prova.xlsx";
data_path = pwd()*"/Data/FinalData/"*data_file;

# Put dates for the start/end of the sample 
start_date = "31/01/1984";
end_date   = "31/12/2024";

# Result Folder name
results_folder = "prova"; # name folder with all the results


# ------------------------------------------------------------------------------
# 3 - Auxiliary Functions
# ------------------------------------------------------------------------------
println("BVAR Set-up > Build auxiliary functions")
eye(x::Int) = Array{Float64,2}(I, x, x);
const j2dt = Dates.julian2datetime;


# ------------------------------------------------------------------------------
# 4 - Load Toolbox Functions 
# ------------------------------------------------------------------------------
println("BVAR Set-up > Load functions in Toolbox")
dir = pwd()*"/Code/Toolbox";
fun = readdir(dir);

for i in fun
    i[end-2:end] == ".jl" ? include(dir*"/"*i) : nothing; 
end


# ------------------------------------------------------------------------------
# 5 - Load Data Section Functions 
# ------------------------------------------------------------------------------
println("BVAR Set-up > Load functions in DataSection")
dir = pwd()*"/Code/DataSection";
fun = readdir(dir);

for i in fun
    i[end-2:end] == ".jl" ? include(dir*"/"*i) : nothing; 
end


# ------------------------------------------------------------------------------
# 5 - Load Data Section Functions 
# ------------------------------------------------------------------------------
println("BVAR Set-up > Load functions in BVAR")
dir = pwd()*"/Code/BVAR";
fun = readdir(dir);
for i in fun
    i[end-2:end] == ".jl" ? include(dir*"/"*i) : nothing; 
end;