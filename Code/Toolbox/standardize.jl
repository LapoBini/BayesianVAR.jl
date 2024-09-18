function standardize(data)

    # --------------------------------------------------------------------------
    # Standardize DataFrame with missing values. Lapo Bini, lbini@ucsd.edu
    # --------------------------------------------------------------------------

    T = size(data,1);
    μ = mean.(skipmissing.(eachcol(data))) |> transpose |> x -> repeat(x, outer = T);
    σ = std.(skipmissing.(eachcol(data)))  |> transpose |> x -> repeat(x, outer = T);;
    X = (data .- μ)./σ;

    return X

end
