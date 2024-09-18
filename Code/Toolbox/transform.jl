function transform(x, flag::Int64; weights::Vector=[], window::Int64=12, index_ma::Vector=[])

    # ------------------------------------------------------------------------------------
    # Applies one of the following transformations:
    #
    # -> flag = 1: window difference
    # -> flag = 2: moving average
    # -> flag = 3: centered moving average
    #
    # Additional arguments can be provided (weights, window, index_ma). It is used within
    # rem_na() function. 
    # ------------------------------------------------------------------------------------

    # initialise
    T = size(x)[1];

    if ~isempty(index_ma)
        output = zeros(size(index_ma)[1]);
    else
        index_ma = collect(1:T);
        output   = zeros(size(x)[1]-window+1);
    end

    if (flag == 3 && window >= size(x)[1])
        error("Error: the dimension of the vector x is lower or equal to the window of the centered MA")
    end

    # window - difference (if window=12, it is an yearly difference)
    if flag == 1
        output = x[1+window:end] - x[1:end-window];
        return output;

    # moving average (simple or centered)
    else
        # if the weights are not specified, take the same weight for all the observations
        if isempty(weights)

            # simple moving average
            if flag == 2
                weights = ones(window)/window;

            # centered moving average
            elseif flag == 3
                weights = ones(convert(Int64, 2*floor(window/2)+1));
                weights = weights./length(weights);
            end
        end

        j = 1;
        for i in index_ma

            # index for calculating the moving average

            # simple moving average
            if flag == 2
                ind_ma = collect(i:-1:i-window+1);

            # centered moving average
            elseif flag == 3
                ind_ma_down = collect(i-1:-1:i-convert(Int64, floor(window/2)));
                ind_ma_up   = collect(i+1:1:i+convert(Int64, floor(window/2)));
                ind_ma      = sort(vcat(ind_ma_down, i, ind_ma_up));
            end

            if (sum(findall(ind_ma.<0)) == 0 && sum(findall(ind_ma.==0)) == 0 && sum(findall(ind_ma.>T)) == 0)
                output[j] = sum(weights[collect(1:size(ind_ma)[1])].*x[ind_ma]);
                j = j + 1;
            end
        end

        return output;
        
    end
end
