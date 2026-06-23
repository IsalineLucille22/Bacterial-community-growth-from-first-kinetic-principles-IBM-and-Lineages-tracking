function Data_daughter = Daughters_Update_fun(Data, val, Lineage_ind)
    Data_temp = Data;
    Data_daughter = Data_temp(Lineage_ind);
    Data_daughter(val == 0) = 0;
end