function cell_ID_plot_daughter = plot_table_fun(cell_ID_plot_daughter, cell_ID_2, ind_Daughters, generation_time_daughter)
temp_div_daughter = reshape(ind_Daughters, [], 2);
ID_D1 = ismember(cell_ID_2, temp_div_daughter(:, 1)); %ID first daughter
ID_D2 = ismember(cell_ID_2, temp_div_daughter(:, 2)); %ID first daughter
cell_ID_plot_daughter(ID_D1) = cell_ID_plot_daughter(ID_D1) + 2.^(-generation_time_daughter(ID_D1));
cell_ID_plot_daughter(ID_D2) = cell_ID_plot_daughter(ID_D2) - 2.^(-generation_time_daughter(ID_D2));
end