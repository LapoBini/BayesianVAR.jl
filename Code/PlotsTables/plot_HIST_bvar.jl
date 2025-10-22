function plot_HIST_bvar(
    IRF::Array{Float64,4},   # Impulse Responses
    HIS::Array{Float64,3},   # Historical Decomposition
    data::DataFrame,         # Vector with data
    H::Int64,                # Forecast horizon IRFs
    name_s::Vector{Any},     # Name of the shock
    p,                       # lag length of the VAR 
    pos_policy,              # Vector with position structural shocks
    base_frq::String,        # Frequency of the BVAR
    results_folder::String,  # where to put results
    transf::Matrix{Bool},    # Transformation to rescale Hist Dec
    ref_dates::Vector{Date}, # In sample dates 
    historical_decomp,       # Variables that I want the hist decomp
    tickers::Vector{Any}
    )

    # --------------------------------------------------------------------------
    # PLOT IRF AND FEVD, Author: Lapo Bini, lbini@ucsd.edu
    # -------------------------------------------------------------------------- 

    # Create Results folder 
    list_dir = readdir("./Results/$results_folder");
    res_path = "./Results/$results_folder/Hist_dec";
    if size(findall(list_dir.==["Hist_dec"]),1) == 0
        mkdir(res_path);
    end

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

    # Confidence interval and alpha for color
    α = [.10, .20, .36]
    c = [0.10, 0.25, 0.85];

    # Loop over the different type of shocks 
    for kk in 1:S

        # Construct median and percentile realization of the shocks  
        e  = mapslices(u -> quantile(u, 0.5), HIS[:,:,kk], dims = 2);
        CU = mapslices(u -> quantile(u, 1.0.-α./2), HIS[:,:,kk], dims = 2);
        CL = mapslices(u -> quantile(u, α./2), HIS[:,:,kk], dims = 2);       

        # Construct ticks for axis 
        data_aux = copy(ref_dates[Int(p[1])+1:end,1]);
        ticks    = DateTime.(unique(year.(data_aux)))[1:10:end];
        tck_n    = Dates.format.(Date.(ticks), "Y");

        # Final djustment 
        ticks[1] = ticks[1] |> lastdayofmonth;
        data_aux = data_aux .|> DateTime; 
        end_tick = data_aux[end];

        # ------------------------------------------------------------------
        # Plot the Realization of the Shock
        # ------------------------------------------------------------------
        Plots.plot(size = (700,500), ytickfontsize  = 13, xtickfontsize  = 13,
                    xguidefontsize = 15, legendfontsize = 13, boxfontsize = 15,
                    framestyle = :box, yguidefontsize = 15, titlefontsize = 18);
        for l in 1:size(CL, 2)
            Plots.plot!(data_aux, CL[:,l], fillrange = CU[:,l], lw = 1, alpha = c[l], 
                        color = "deepskyblue1", label = "")
        end
        hline!([0], color = "black", lw = 1, label = nothing)
        Plots.plot!(data_aux, e[:], lw = 3, color = "black", label = "")
        Plots.plot!(ylabel = "", title = name_s[kk]*" - Historical Decomposition",
                    left_margin = 1Plots.mm, right_margin = 1Plots.mm,
                    bottom_margin = 1Plots.mm, top_margin = 1Plots.mm,
                    xticks = (ticks,tck_n), xlims = (ticks[1], end_tick))
        Plots.savefig(res_path*"/$(name_s[kk]).pdf")

        # ------------------------------------------------------------------
        # Plot the contribution to the instrumented variables
        # ------------------------------------------------------------------
        # Initialize excel file 
        intro     = [["";""] [""; ""]];
        est_HIST  = DataFrame(intro, Symbol.([" ",""]))
        res_excel = res_path*"/$(replace(name_s[1], " "=>""))_HIST_.xlsx";

        # Allocate empty spreadsheet
        XLSX.openxlsx(res_excel, mode = "w") do file
                    
            # Save FEVD 
            XLSX.rename!(file[1], "EMPTY")
            XLSX.writetable!(file[1], est_HIST)

        end

        # Find position shock of interest
        pos = findall(pos_policy .== kk)[1]

        for hd in 1:length(historical_decomp)
            ii  = findall(tickers .== historical_decomp[hd])[1];
            Ω   = mapslices(u -> quantile(u, 0.5), IRF[:,ii,:,pos], dims = 2);
            cu  = mapslices(u->quantile(u, 1.0.-α./2), IRF[:,ii,:,pos], dims = 2);
            cl  = mapslices(u->quantile(u, α./2), IRF[:,ii,:,pos], dims = 2);
            
            # Obtain BC component by rescaling by Real GDP IRFs
            (transf[ii,1] == 1) .& (transf[ii,2] == 0) .& (transf[ii,3] == 0) ? scale = 100 : scale = 1;
            T    = size(e, 1);
            irfᵢ = kron(e', Ω.*scale); 
            BC   = zeros(T);

            for t in 1:T
                aux = [irfᵢ[t-(i-1),i] for i in 1:t]
                BC[t] = sum(aux)
            end

            # Ontain confidence interval (based on the median response)
            AUX = HIS[:,:,kk]
            BCα = zeros(T,size(AUX,2))

            for i in 1:size(AUX,2)
                irfᵢ = kron(AUX[:,i]', Ω.*scale); 
                for t in 1:T
                    aux = [irfᵢ[t-(i-1),i] for i in 1:t]
                    BCα[t,i] = sum(aux)
                end
            end

            # Confidence interval historical decomposition
            cu = mapslices(u -> quantile(u, 1.0.-α./2), BCα, dims = 2);
            cl = mapslices(u -> quantile(u, α./2), BCα, dims = 2);

            # Plot
            name_p = "$(replace(name_s[1], " "=>""))_$(historical_decomp[hd])"
            Plots.plot(size = (700,500), ytickfontsize  = 13, xtickfontsize  = 13,
                        xguidefontsize = 15, legendfontsize = 13, boxfontsize = 15,
                        framestyle = :box, yguidefontsize = 15, titlefontsize = 18);
            for l in 1:length(α)
                Plots.plot!(data_aux, cl[:,l], fillrange = cu[:,l], lw = 1, alpha = c[l], 
                            color = "deepskyblue1", label = "")
            end
            hline!([0], color = "black", lw = 1, label = nothing)
            Plots.plot!(data_aux, BC[:], lw = 3, color = "black", label = "")
            Plots.plot!(ylabel = "", title = names(data)[ii+1],
                        left_margin = 1Plots.mm, right_margin = 1Plots.mm,
                        bottom_margin = 1Plots.mm, top_margin = 1Plots.mm,
                        xticks = (ticks,tck_n), xlims = (ticks[1], end_tick))
            Plots.savefig(res_path*"/$(name_p).pdf")

            # Save file with cycle component 
            col_names_excel = "Hist" .* [""; string.(([1.0.-α./2; α./2] .|> u->round(u, digits = 2)))]
            aux_df = DataFrame([BC cl cu], Symbol.(col_names_excel))

            # Open excel file and add extra spreadsheet
            XLSX.openxlsx(res_excel, mode = "rw") do file

                # Add extra sheet 
                sheet = XLSX.addsheet!(file, historical_decomp[hd]);
                XLSX.writetable!(sheet, aux_df)
            end
        end
    end
end
