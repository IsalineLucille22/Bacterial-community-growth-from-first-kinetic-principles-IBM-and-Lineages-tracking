function Data = Correction_Cell_Sizes(Data, Length_to_corr, mean_val, fact_min, fact_max)
Data.(Length_to_corr) = max(min(Data.(Length_to_corr), fact_max*mean_val), fact_min*mean_val);
end