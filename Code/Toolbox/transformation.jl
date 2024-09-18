function transformation(
    A,                  # Variables to transform
    t,                  # flags for transformation
    freq::String;       # frequency variable 
    arithmetic_g = true # arithmetic or logarithmic transf
    )

# ------------------------------------------------------------------------------
# Author: Lapo Bini, lapo.bini@barcelonagse.eu
# ------------------------------------------------------------------------------

    nvars = size(t,1)
    X     = copy(A);

    for vars in 1:nvars

        # If you want arithmetic YoY growth
        if (sum(t[vars,:]) .== 3) && arithmetic_g

            # Frequency of the variables
            if freq == "d"
                l = 264
            elseif freq =="w"
                l = 52
            elseif freq == "m"
                l = 12
            elseif freq == "q"
                l = 4
            end

            xxx = A[l+1:end,vars]
            for i in 1:length(xxx)
                xxx[i] = ((xxx[i] - A[i,vars])/A[i,vars])*100
            end
            X[1+l:end, vars] =  xxx;
            X[1:l,vars]     .= missing;

        else
            # log
            if t[vars,1] == true
                X[:,vars] = log.(X[:,vars])
            end

            # first difference
            if t[vars,2] == true
                X[:,vars] = [missing; X[2:end,vars]-X[1:end-1,vars]]
            end

            # MoM, QoQ or WoW growth rate depending on the frequency of A 
            if (t[vars,1] == true) && (t[vars,2] == true)
                X[:,vars] = 100*X[:,vars];
            end

            # Filters: if you want logarithmic YoY or diff YoY growth
            if (t[vars,2] == true) && (t[vars,3] == true)
                if freq == "d"
                    l = 264
                elseif freq =="w"
                    l = 52
                elseif freq == "m"
                    l = 12
                elseif freq == "q"
                    l = 4
                end

                # Our case Flag = 2 because we alredy took log and diff.
                w1        = ones(l);
                X[:,vars] = use_filter(X[:,vars], 2, weights=w1, window=l);
            end
        end
    end

    return X

end
