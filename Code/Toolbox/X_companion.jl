function X_companion(
    y::Array{Float64,2},
    p::Int64;
    intercept_c = true
    )

    t,k = size(y);
    idx = collect(1:k:k*p+1)
    X   = zeros(k*p, t-p)

    @inbounds for j = 1:p
        X[idx[j]:idx[j+1]-1,:] = y[p+1-j:t-j,:]'
    end

    if intercept_c
        X = [X; ones(1, t-p)]
    end

    # Dimension of the output is (k⋅p + 1)ₓ(t-p) if we include the intercept
    return X

end
