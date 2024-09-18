function use_filter(x, flag::Int64; weights::Vector=[], window::Int64=[], index_ma::Vector=[])

    #=
    --------------------------------------------------------------------------------------------------------------------
    Description:   Applies one of the following filters to x and returns the outcome:
                   a) flag = 1: difference (window)
                   b) flag = 2: moving average (window)
                   c) flag = 3: centered moving average (window)

                   Additional arguments can be used (weights, window, index_ma)
    --------------------------------------------------------------------------------------------------------------------
    =#
    # Dimensions
    m = size(x)[1];

    # Error management
    if (flag == 3) && (window >= m)
        error("Error: the dimension of the vector x is lower or equal to the window of the centered MA")
    end

    if isempty(index_ma)
        index_ma = [1:m;];#collect(1:m);
    end

    if flag == 1
        output = zeros(size(x)[1]);

    elseif (flag == 2) || (flag == 3)
        output = zeros(size(index_ma)[1]);
    end

    if flag == 3
        window_plus_1 = convert(Int64, 2*floor(window/2)+1);
        half_window   = convert(Int64, floor(window/2));
    end

    if isempty(window)
        window = convert(Int64, sum(weights));
    end

    # ----- Execution -----

    # Difference (e.g., window=12 -> yearly difference)
    if flag == 1
        output[1:window]     .= missing;
        output[window+1:end] = x[1+window:end]-x[1:end-window];
        return output;

    # Moving average (simple or centered)
    else

        # If the weights are not specified, take the same weight for all the observations
        if isempty(weights)

            # Simple moving average
            if flag == 2
                weights = ones(window)/window;

            # Centered moving average
            elseif flag == 3
                weights = ones(window_plus_1)./window_plus_1;
            end
        end

        # Iterate index_ma
        if flag == 2
            j = window;
            if window > 1
                output[1:window-1] .= missing;
            end

        elseif flag == 3
            j = half_window+1;
            if half_window > 1
                output[1:half_window]         .= missing;
                output[end-half_window+1:end] .= missing;
            end
        end

        for i in index_ma

            # Index: simple moving average
            if flag == 2
                ind_ma = collect(i:-1:i-window+1);

            # Index: centered moving average
            elseif flag == 3
                ind_ma_down = collect(i-1:-1:i-half_window);
                ind_ma_up   = collect(i+1:1:i+half_window);
                ind_ma      = sort(vcat(ind_ma_down, i, ind_ma_up));
            end

            if (sum(findall(ind_ma.<0)) == 0) && (sum(findall(ind_ma.==0)) == 0) && (sum(findall(ind_ma.>m)) == 0)
                output[j] = sum(weights[collect(1:size(ind_ma)[1])].*x[ind_ma]);
                j += 1;
            end
        end

        return output;
    end
end
