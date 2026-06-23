function Data = Remove_Cells(Data, Length_to_corr, mean_val, fact_min, fact_max)
ind_remove = Data.(Length_to_corr) <= fact_min*mean_val;
Data(ind_remove, :) = [];
ind_remove = Data.(Length_to_corr) >= fact_max*mean_val;
Data(ind_remove, :) = [];
end