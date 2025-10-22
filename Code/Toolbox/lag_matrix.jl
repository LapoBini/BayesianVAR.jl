function lag_matrix(X::Matrix, p::Int)

    # Dimensions of input matrix
    T, N = size(X)
    
    # Total number of columns in the output matrix: N columns per lag, 
    # p+1 sets (contemporaneous + p lags)
    num_cols = N * (p + 1)
    
    # Initialize output matrix
    lagged_matrix = zeros(T-p, num_cols)
    
    # Fill the matrix  create lagged value shock 
    for i in 0:p
        col_start = (i * N) + 1
        col_end   = (i + 1) * N
        lagged_matrix[:,col_start:col_end] = X[p+1-i:end-i,:]
    end
    
    return lagged_matrix

end
