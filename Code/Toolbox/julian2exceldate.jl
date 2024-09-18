function julian2exceldate(julian_date)
    
    excel_date = julian_date.-Dates.datetime2julian(DateTime(1899,12,30));

    return excel_date
end
