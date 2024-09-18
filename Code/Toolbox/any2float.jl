function any2float(A)

    # Convert any datatype with missing to Array{Float64,2}
    ncols = size(A,2)
    nrows = size(A,1)

    # Preallocate 
    B = Array{Float64}(zeros(size(A)))

    for rows in 1:nrows
        for cols in 1:ncols
            if ismissing(A[rows,cols])
                B[rows,cols] = NaN
            else
                B[rows,cols] = A[rows,cols]
            end
        end
    end

    return B
end
