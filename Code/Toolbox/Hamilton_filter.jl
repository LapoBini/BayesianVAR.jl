function Hamilton_filter(
    y::Vector{Any};
    h = 24,
    p = 12
    )

    # Create lagged and contemporaneous values 
    X   = lag_matrix(y[:,:], p) |> Matrix{Float64};
    yₜ₊ₕ = X[h+1:end,1];
    T   = length(yₜ₊ₕ);
    x   = [ones(T,1) X[1:end-h,2:end]];

    # Compute parameter 
    α = (x' * x) \ (x' * yₜ₊ₕ);
    τ = x * α;
    c = yₜ₊ₕ - τ;

    return [zeros(p+h,2) .* NaN; τ c]

end