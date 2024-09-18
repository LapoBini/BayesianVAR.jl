
function readdata_haver(
    # Partial path to the excel file
    data_path::String,
    # Dates for in-sample and out-of-sample estimation
    start_date::String,
    end_date::String;
    save_transformed = false
    )

    # --------------------------------------------------------------------------
    # LOAD DATASET, Author: Lapo Bini, lbini@ucsd.edu
    # -------------------------------------------------------------------------- 
    # REMEMBER :
    # 1 - It loads automatically all the different frequencies (and the use 
    #     interpolation). The important thing is to create a unique legend and 
    #     sign worksheets where the Highest frequency variables always match the 
    #     order in the Data secion. 
    # 2 - Call the worksheets with data "Data_" + frequency ("d", "w", "m", "q").
    # 3 - If using sign identification, 1 for positive response, -1 for negative 
    #     responses. The first column of "Sign" worksheet must have the position 
    #     of the policy variable affected by the shock (for example one type of 
    #     sign restrictio, then you put 1 on the variable affected by the shock,
    #     if two shocks you must specify 1 for the first shock, 2 for the second)
    # 4 - If using instrumental variable, use a worksheet "IV" to upload the 
    #     instrument. Remember also to add the dates of the shocks.  
    # --------------------------------------------------------------------------

    # --------------------------------------------------------------------------
    # Load File and Save Auxiliary Stuff 
    # --------------------------------------------------------------------------
    # Useful for aggregation 
    dict_agg = Dict("d" => 1, "w" => 2, "m" => 3, "q" => 4, "y" => 5);

    # Load excel file 
    database  = XLSX.readxlsx(data_path);
    sheetname = XLSX.sheetnames(database);

    # Create all auxiliary variables to transofrm data and other stuff 
    # Included series 
    legend = database["Legend"][:][2:end,:];
    incl   = legend[:,end] |> Array{Bool,1};

    # For variables name: put the name that you wnat in the second column of the
    # spreadshit Legend
    var_names = legend[incl,2];
    prior     = legend[incl,end-1] |> Array{Float64,1};
    transf    = legend[:, 3:5] |> Array{Bool,2};
    tickers   = legend[incl,1];
    for i in axes(tickers,1)
        tickers[i] = filter(x -> !isspace(x), tickers[i])
    end;

    # Transform data 
    pos = [sheetname[i][1:2] for i in 1:length(sheetname)] |> x -> findall(x .== "Da");
    frq = [sheetname[i][end] for i in pos] .|> string;


    # --------------------------------------------------------------------------
    # Transformation Series
    # --------------------------------------------------------------------------
    dic_data = Dict();
    dic_date = Dict();
    count    = 0;

    for i in 1:length(pos)

        # Select data 
        aux     = database[sheetname[pos[i]]][:];
        X       = aux[19:end,3:end];
        N       = size(X,2)

        # Transform and save data 
        idx_aux  = count+1:count+N; # Select position variable given the frequency 
        incl_aux = incl[idx_aux];   # including boolean for the frequency 
        X_transf = transformation(X[:,incl_aux], transf[idx_aux,:][incl_aux,:], frq[i]);

        # Save insample dates
        ref_dates = aux[19:end, 2] .|> Date;

        # Save transformed data
        dic_data[frq[i]] = X_transf;
        dic_date[frq[i]] = ref_dates;

        # Update count 
        count += N;

    end


    # --------------------------------------------------------------------------
    # Aggregation to High Frequency
    # --------------------------------------------------------------------------
    # Find highest Frequency
    agg      = [dict_agg[i] for i in frq];
    base_frq = argmin(agg);

    # Select insample dates and frequency that we need to interpolate and put them in 
    # order from the higher to the lower frequency
    order_frq = ["d", "w", "m", "q", "y"];
    frq_aux   = frq[1:end .!= base_frq];
    frq_aux   = [frq_aux[findall(frq_aux .== order_frq[i])] for i in 1:length(order_frq)] |> x->vcat(x...);

    # Select the highest frequency, that would be the baseline frequency if the model 
    base_frq  = frq[base_frq];
    ref_dates = dic_date[base_frq];
    y_high    = dic_data[base_frq];

    # Create final Array and allocate highest frequency variables 
    Y = Array{Any}(missing, length(ref_dates), sum(incl));
    c = size(y_high,2); 
    Y[:,1:c] = y_high;

    for i in 1:length(frq_aux)

        # Load data and date
        y_aux = dic_data[frq_aux[i]];
        d_aux = dic_date[frq_aux[i]];

        # Allign date if lower frequency starts earlier 
        start_low = findall(d_aux .>= ref_dates[1])[1]

        # High frequency counterpart with missings
        y_aux_high, _ = low2highfrequency(y_aux[start_low:end,:], ref_dates, d_aux[start_low:end]);

        # Interpolate missings (considering leading and closing missings)
        y_aux_int, _, ind_lead_end = rem_na(y_aux_high, option=0, k = 3)
        y_aux_high[ind_lead_end,:] = y_aux_int;

        # Allocate
        c1 = size(y_aux_high,2);
        Y[:,c+1:c+c1] = y_aux_high;
        c += c1
        
    end

    # Take balanced sample 
    data, ind_mis, ind_balance = balanced_sample(Y);
    ref_dates   = ref_dates[ind_balance];

    # Select the desired priod
    start_d     = Date(start_date, "dd/mm/yyyy")
    end_d       = Date(end_date, "dd/mm/yyyy")
    sample_bool = (ref_dates .>= start_d) .& (ref_dates .<= end_d)
    ref_dates   = ref_dates[sample_bool];
    data        = data[sample_bool,:];
    transf_data = DataFrame([ref_dates data], Symbol.(["Dates"; var_names]));

    # --------------------------------------------------------------------------
    # Identification - Sign Restriction or IV?
    # --------------------------------------------------------------------------
    # if 1 positive response, if -1 negative response, if missing unrestricted. 
    # Column name must be the name of the shock  
    sign_s     = [];
    name_s     = [];
    pos_policy = [];
    if "Sign" in sheetname
        sign_s     = database["Sign"][:][2:end,3:end][incl,:] |> any2float;
        name_s     = database["Sign"][:][1,3:end];
        pos_policy = database["Sign"][:][2:end,1][incl] |> any2float;
    end

    # External instrument. Remember the column name must be the ticker of the 
    # series as it is in the Legend spreadsheet. First column in the variable
    # instrumeny will be the date of the instrument
    instrument   = [];
    instrumented = [];
    if "IV" in sheetname
        iv_aux     = database["IV"][:];
        instrument = DataFrame(iv_aux[2:end,1:end], Symbol.(iv_aux[1,1:end]));
    end


    # --------------------------------------------------------------------------
    # Save Final Dataset and Return 
    # --------------------------------------------------------------------------
    if save_transformed == true;
        res   = readdir(pwd()*"/Data");
        res_f = "FinalData";
        res_f in res ? nothing : mkdir(pwd()*"/Data/"*res_f);

        XLSX.openxlsx(pwd()*"/Data/"*res_f*"/transformed_data.xlsx", mode="w") do file
            
            # Save final data (first spreedshet has been already created)
            XLSX.rename!(file[1], "Data")
            XLSX.writetable!(file[1], transf_data)

        end
    end

    return transf_data, ref_dates, tickers, prior, sign_s, name_s, 
           instrument, pos_policy, base_frq

end