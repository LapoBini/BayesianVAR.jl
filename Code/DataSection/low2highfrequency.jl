function low2highfrequency(a, ref_dates, ref_freq)
    
    # --------------------------------------------------------------------------
    # Author: Lapo Bini, lbini@ucsd.edu
    # --------------------------------------------------------------------------
    
    # Ref_freq is for the low frequency, ref_dates are the in-sample dates 
    # for the high frequency variables 
    N_freq        = size(a,2);
    b             = Array{Any}(missing, length(ref_dates), N_freq);
    idx_freq = Array{Int64}(zeros(length(ref_freq)));

    for i = 1:length(ref_freq)
        idx_freq[i] = maximum(findall(ref_dates .<= ref_freq[i]));
    end

    idx_freq = unique(idx_freq);
    ind_seq  = collect(1:1:length(ref_freq));

    b[idx_freq, :] = a[ind_seq, :];

    return b, idx_freq

end
