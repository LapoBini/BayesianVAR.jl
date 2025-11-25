function plot_TREND_bvar(
    TR::Array{Float64,3},    # trend
    data::DataFrame,         # Vector with data
    p,                       # lag length of the VAR 
    results_folder::String,  # where to put results
    ref_dates::Vector{Date}, # In sample dates 
    tickers::Vector{Any}
    )

    # --------------------------------------------------------------------------
    # Plot Trend of the VAR, Author: Lapo Bini, lbini@ucsd.edu
    # -------------------------------------------------------------------------- 

    # Create Results folder 
    list_dir = readdir("./Results/$results_folder");
    res_path = "./Results/$results_folder/trend";
    if size(findall(list_dir.==["trend"]),1) == 0
        mkdir(res_path);
    end

    # Name variables 
    var_names = names(data)[2:end];
    k         = length(var_names)

    # Confidence interval and alpha for color
    α  = [.10, .20, .36]
    c  = [0.10, 0.25, 0.85];
    cc = [RGB(0, 0.4470, 0.7410); RGB(0.8500, 0.3250, 0.0980)];

    # Initialize excel file 
    intro     = [["";""] [""; ""]];
    est_TR    = DataFrame(intro, Symbol.([" ",""]))
    res_excel = res_path*"/trends.xlsx";

    # Allocate empty spreadsheet
    XLSX.openxlsx(res_excel, mode = "w") do file
                    
        # Save trend
        XLSX.rename!(file[1], "EMPTY")
        XLSX.writetable!(file[1], est_TR)
    end

    # Loop over the different type of shocks 
    for kk in 1:k

        # Construct median and percentile realization of the shocks  
        e  = mapslices(u -> quantile(u, 0.5), TR[kk,:,:], dims = 2);
        CU = mapslices(u -> quantile(u, 1.0.-α./2), TR[kk,:,:], dims = 2);
        CL = mapslices(u -> quantile(u, α./2), TR[kk,:,:], dims = 2);       

        # Construct ticks for axis 
        data_aux = copy(ref_dates[Int(p[1])+1:end,1]);
        ticks    = DateTime.(unique(year.(data_aux)))[1:10:end];
        tck_n    = Dates.format.(Date.(ticks), "Y");

        # Final djustment 
        ticks[1] = ticks[1] |> lastdayofmonth;
        data_aux = data_aux .|> DateTime; 
        end_tick = data_aux[end];

        # ------------------------------------------------------------------
        # Plot Trend and Actual Series
        # ------------------------------------------------------------------
        Plots.plot(size = (700,500), ytickfontsize  = 13, xtickfontsize  = 13,
                    xguidefontsize = 15, legendfontsize = 13, boxfontsize = 15,
                    framestyle = :box, yguidefontsize = 15, titlefontsize = 18);
        for l in 1:size(CL, 2)
            Plots.plot!(data_aux, CL[:,l], fillrange = CU[:,l], lw = 1, alpha = c[l], 
                        color = "deepskyblue1", label = "")
        end
        Plots.plot!(data_aux, e[:], lw = 3, color = "black", label = "Trend")
        Plots.plot!(data_aux, data[Int(p[1])+1:end, 1+kk], lw = 3, color = cc[2],
                    label = "Actual")
        Plots.plot!(ylabel = "", title = var_names[kk],
                    left_margin = 1Plots.mm, right_margin = 1Plots.mm,
                    bottom_margin = 1Plots.mm, top_margin = 1Plots.mm,
                    xticks = (ticks,tck_n), xlims = (ticks[1], end_tick))
        names_p = filter(x -> !isspace(x), var_names[kk])
        Plots.savefig(res_path*"/$names_p.pdf")

        # ------------------------------------------------------------------
        # Save Trend in Excel File 
        # ------------------------------------------------------------------
        # Save file with cycle component 
        col_names_excel = "TR" .* [""; string.(([1.0.-α./2; α./2] .|> u->round(u, digits = 2)))]
        aux_df = DataFrame([data_aux e CL CU], Symbol.(["date"; col_names_excel]))

        # Open excel file and add extra spreadsheet
        XLSX.openxlsx(res_excel, mode = "rw") do file

            # Add extra sheet 
            sheet = XLSX.addsheet!(file, tickers[kk]);
            XLSX.writetable!(sheet, aux_df)
        end
    end
end
