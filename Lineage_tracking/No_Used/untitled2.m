daugt_test = [Table_fin_tot.Daughter_ID_1; Table_fin_tot.Daughter_ID_2];
daugt_test = daugt_test(daugt_test ~= 0);
unique_vals_Daughters = unique(daugt_test);
more_two_daughters = histcounts(daugt_test, [unique_vals_Daughters; max(unique_vals_Daughters) + 1]);