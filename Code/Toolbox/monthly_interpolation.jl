function monthly_interpolation(X, k)

    t, n  = size(X);
    X_res = Array{Union{Missing, Float64},2}(missing, t,n);

    for i in 1:n

        non_missing  = findall(.~ismissing.(X[:,i]));
        length_miss  = non_missing[2:end] - non_missing[1:end-1];
        non_missing2 = [non_missing[findall(length_miss .<= k)]; non_missing[end]];
        interpolate  = X[non_missing2[1]:non_missing2[end],i:i]

        x_interpolated,_,_ = rem_na(interpolate, option=0, k = 3, restricted = false)
        X_res[non_missing2[1]:non_missing2[end],i] = x_interpolated;

        X_res[non_missing2,i] = X[non_missing2,i]
    end

    return X_res

end
