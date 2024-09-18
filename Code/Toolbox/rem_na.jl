function rem_na(y; option::Int64=1, k::Int64=3)

    # ----------------------------------------------------------------------------------------
    # Handle the NAs in y using one out of two methods:
    #
    # option == 0 -> replace missing values after removing leading and closing zeros
    # option == 1 -> only remove rows with leading and closing zeros (default)
    #
    # 2*k+1 is the window for the simple moving average required if option == 1 
    # ----------------------------------------------------------------------------------------

    T, n   = size(y);
    ind_na = ismissing.(y);        

    if option == 0
        ind_rem = sum(ind_na, dims=2) .> n*0.8;
    elseif option == 1
        ind_rem = sum(ind_na, dims = 2) .== n;
    end

    ind_na_lead = cumsum(ind_rem, dims=1) .== collect(1:T);
    ind_na_end  = cumsum(ind_rem[end:-1:1], dims=1) .== collect(1:T);
    ind_na_end  = ind_na_end[end:-1:1];

    
    ind_not_lead_end_na = findall(sum(ind_na_lead, dims = 2)[:,1] + sum(ind_na_end, dims = 2)[:,1] .!= 1);
    y                   = y[ind_not_lead_end_na, :]; # remove leading and closing NAs
    ind_na              = ismissing.(y); # update ind_na

    if option == 0
        for i in 1:n

            # replace missing values after removing leading and closing NAs
            x             = y[:, i];
            ind_is_not_na = findall(.~(ismissing.(x)));

            t1            = ind_is_not_na[1];
            t2            = ind_is_not_na[end];
            spl           = Spline1D(ind_is_not_na, x[ind_is_not_na]);
            x[t1:t2]      = spl(t1:t2);

            # This removes closing NAs for y-ith, whenever you don't have NAs for the other variables
            ind_na_x = ismissing.(x);

            if sum(ind_na_x) > 0
                x[ind_na_x] .= median(x[.~ind_na_x]);

                # replace all entries at i, where i >= window (2*k+1)
                ind_nan_x_ma = findall(ind_na_x);

                if ~isempty(ind_nan_x_ma)

                    x_MA = transform((vcat(x[1]*ones(k), x, x[end]*ones(k))), 2,
                                        weights=ones(2*k+1)./(2*k+1), window=2*k+1, index_ma=ind_nan_x_ma.+2*k);

                    x[ind_nan_x_ma] .= x_MA;
                end
            end

            y[:, i] = x;
        end
    end

    return y, ind_na, ind_not_lead_end_na;
end
