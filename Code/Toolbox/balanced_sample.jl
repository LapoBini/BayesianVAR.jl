function balanced_sample(data)

    # --------------------------------------------------------------------------
    # Remove Leading/Closing missing,  Author: Lapo Bini, lbini@ucsd.edu
    # --------------------------------------------------------------------------

    # Find missing values 
    T       = size(data, 1);
    ind_mis = ismissing.(data) |> Array{Bool};
    aux     = sum(ind_mis, dims = 2) .> 0;

    # ind_na_lead: As soon as no missing, cumsum(1:t) ≂̸ t 
    # ind_na_end: same but for closing missing 
    ind_na_lead = cumsum(aux, dims=1) .== collect(1:T);
    ind_na_end  = (cumsum(aux[end:-1:1], dims=1) .== collect(1:T))[end:-1:1];
    ind_balance = findall(sum(ind_na_lead, dims=2)[:] + sum(ind_na_end, dims = 2)[:] .!= 1);

    return data[ind_balance,:], ind_mis, ind_balance

end
