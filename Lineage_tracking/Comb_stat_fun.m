function cellcounts = Comb_stat_fun(cellcounts, new_table, table_temp, sub_table_temp, stat_tp, name_table, i, ind_1, val_temp, temp_stat)
    t_1 = table_temp.Timepoint(1);
    new_table(:, (t_1 + 1):stat_tp) = sub_table_temp(:, 1:(stat_tp - t_1));
    new_table(val_temp, (t_1 + 1):end) = [sub_table_temp(val_temp, 1:(stat_tp - t_1)) temp_stat(ind_1, 2:end)];
    cellcounts.lineage(i).(name_table) = new_table;
end