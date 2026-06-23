function Data = Arrage_Timepoint_fun(Data)
unique_time = unique(Data.Timepoint);
new_unique_time_point = 0:length(unique_time) - 1;
Data_temp = Data.Timepoint;
for i = 1:length(unique_time)
    % ind_temp = find(Data.Timepoint == unique_time(i));
    Data_temp(Data.Timepoint == unique_time(i)) = new_unique_time_point(i);
end
Data.Timepoint = Data_temp;
end