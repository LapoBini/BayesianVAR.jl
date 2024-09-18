function excel2datetime(A)

    ncols = size(A,2)
    nrows = size(A,1)

    for row in 1:nrows
        for col in 1:ncols
            if isnan(A[row,col]) == false
                A[row,col] = A[row,col] + Dates.datetime2julian(DateTime(1899,12,30))
            end
        end
    end
    
    return A
end
