function quarterly2monthly(
    Xq,      # vector or matrix to transform to monthly frequency 
    n::Int64 # length of the monthly vector that we want 
    )
    
    Xqm = kron(Xq, [missing; missing; 1]);
    if ndims(Xqm) == 1
        Xqm = [Xqm; missing.*ones(n-size(Xqm,1))];
    else
        Xqm = [Xqm; missing.*ones(n-size(Xqm,1), size(Xqm,2))];
    end

    return Xqm;

end
