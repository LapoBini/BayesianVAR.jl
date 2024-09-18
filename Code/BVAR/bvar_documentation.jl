function bvar_documentation(
    IRF::Array{Float64,4}, # Impulse Responses
    data::DataFrame,       # Vector with data
    H::Int64,              # Forecast horizon IRFs
    name_s,                # Name of the shock
    pos_policy,            # Vector with position structural shocks
    base_frq::String,      # Frequency of the BVAR
    res_path::String       # where to put results
    )

    # --------------------------------------------------------------------------
    # PLOT IRF AND FEVD, Author: Lapo Bini, lbini@ucsd.edu
    # -------------------------------------------------------------------------- 

    # Name variables 
    var_names = names(data)[2:end];
    k         = length(var_names)
    Hᵢ        = size(IRF,1);

    # Number of structural shocks considered
    if typeof(pos_policy) == Int64 # Case for IV shock
        pos_S = [pos_policy];
        S     = length(pos_S)
    else # Case for Sign restriction
        pos_S = findall(.~isnan.(pos_policy))
        S     = length(pos_S);
    end

    # Auxiliary variables for plotting 
    x_ax  = collect(0:1:H);
    if base_frq == "m"
        x_label = "Months"
        ticks   = [0; collect(12:12:H)]
    else
        x_label = "Quarters"
        ticks   = [0; collect(4:4:H)]
    end

    # Confidence interval and alpha for color
    α     = [.10, .20, .36]
    c     = [0.10, 0.25, 0.85];

    
    # ----------------------------------------------------------------------
    # LOOP OVER TYPE OF SHOCKS 
    # ----------------------------------------------------------------------
    for i in 1:S

        # Create directory a precise shock
        title_i   = name_s[i];
        title_dir = replace(title_i,  " " => "_")
        list_dir  = readdir(res_path);
        if size(findall(list_dir .== title_dir),1) == 0
            mkdir(res_path*"/"*title_dir );
        end


        # ------------------------------------------------------------------
        # 1 - Plot Impulse Response Functions 
        # ------------------------------------------------------------------
        # Allocate confidence interval
        CU = Array{Float64,3}(undef, Hᵢ, length(α), k);
        CL = Array{Float64,3}(undef, Hᵢ, length(α), k);
        IR = Array{Float64,2}(undef, Hᵢ, k);

        # Compute the quantiles from the simulated IRFs
        typeof(pos_policy) == Int64 ? pos_aux = [1] : pos_aux = pos_S
        for j = 1:k
            CU[:,:,j] = mapslices(u->quantile(u, 1.0.-α./2), IRF[:,j,:,pos_aux[i]], dims = 2);
            CL[:,:,j] = mapslices(u->quantile(u, α./2), IRF[:,j,:,pos_aux[i]], dims = 2);
            IR[:,j]   = mapslices(u->quantile(u, 0.5), IRF[:,j,:,pos_aux[i]], dims = 2);
        end

        println("BVAR > Structural Analysis > Plot results > IRFs")
        b₀ᵢ = IR[1,pos_S[i]];
        for j in 1:k 
            name_p = filter(x -> !isspace(x), var_names[j]);
            Plots.plot(size = (700,500), ytickfontsize  = 15, xtickfontsize  = 15,
                        xguidefontsize = 15, legendfontsize = 13, boxfontsize = 15,
                        framestyle = :box, yguidefontsize = 15, titlefontsize = 25);
            for l in 1:size(CL[:,:,1], 2)
                Plots.plot!(x_ax, CL[:,l,j]./b₀ᵢ, fillrange = CU[:,l,j]./b₀ᵢ,
                            lw = 1, alpha = c[l], color = "deepskyblue1", xticks = ticks,
                            label = "")
            end
            Plots.plot!(x_ax, IR[:,j]./b₀ᵢ, lw = 3, color = "black", xticks = ticks,
                        label = "")
            hline!([0], color = "black", lw = 1, label = nothing)
            Plots.plot!(xlabel = "Months", ylabel = "", title = var_names[j],
                        left_margin = 1Plots.mm, right_margin = 3Plots.mm,
                        bottom_margin = 1Plots.mm, top_margin = 7Plots.mm,
                        xlims = (0,H))
            Plots.savefig(res_path*"/"*title_dir*"/"*string(name_p)*".pdf")
        end
    end
end
