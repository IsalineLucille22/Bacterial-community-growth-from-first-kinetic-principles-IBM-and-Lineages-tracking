function [Data_kept_init, translated_vector] = AdditionNoise(Data_kept_init)
time_point = unique(Data_kept_init.Timepoint);
translated_vector = zeros(1, length(time_point));
diff_x = 0; diff_y = 0;
for i = 2:length(time_point)
    ind_temp = find(Data_kept_init.Timepoint == time_point(i));
    diff_x = diff_x + normrnd(0, 0.1); diff_y = diff_y + normrnd(0, 0.1); 
    Data_kept_init.Centroid_x(ind_temp) = Data_kept_init.Centroid_x(ind_temp) + diff_x; %Addition of a translation movement 
    Data_kept_init.Centroid_y(ind_temp) = Data_kept_init.Centroid_y(ind_temp) + diff_y;
    translated_vector(i) = sqrt(diff_x^2 + diff_y^2);
end
end