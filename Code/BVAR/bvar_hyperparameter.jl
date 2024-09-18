@with_kw mutable struct Hyperparameter
    λ::Array{Float64,1}
    τ::Array{Float64,1}
    γ::Array{Float64,1}
    ε::Float64
    p::Array{Int64,1}
    H::Int64
    reps::Int64
    burnin::Int64
    max_try::Int64
    update::Int64
end
