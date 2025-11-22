function VAR_trend(
    A::AbstractMatrix, 
    x0::AbstractVector, 
    t::Int
    )

    n   = size(A, 2)            
    lag = Int((size(A, 1)-1)/n)    

    # Companion matrix F (12×12 here)
    F = [A[1:end-1,:]'; eye(n*(lags-1)) zeros(n*(lags-1), n)];

    # Constant in companion form (zero by default)
    C = zeros(n*lag);  
    C[1:n] = A[end,:]

    τₜ = zeros(eltype(A), n, t)
    s = copy(x0)
    for t in 1:t
        s = F * s + C
        τₜ[:, t] = s[1:n]
    end

    return τₜ
    
end
