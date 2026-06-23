function  Data_kept = Length_cells_mod(Data_kept, mean_L, fact_Thres)
 Length_cell = Data_kept.axis_major_length;
 ind_mod = find(Length_cell > fact_Thres*mean_L);
 Length_cells_to_mod = Length_cell(ind_mod);
 Pos_X = [Data_kept.Centroid_x(ind_mod), Data_kept.Centroid_y(ind_mod)]';
 angle_cell = Data_kept.orientation(ind_mod);
 New_ID = (max(Data_kept.Mask_nb) + 1):(max(Data_kept.Mask_nb) + length(ind_mod));
 for i = 1:length(ind_mod)
    x_center = Pos_X(1, i);
    y_center = Pos_X(2, i);
    L_seg = Length_cells_to_mod(i)/4;
    theta = angle_cell(i);
    L_seg_cos = L_seg*cos(theta); L_seg_sin = L_seg*sin(theta);
    Pos_X_1 = [L_seg_cos + x_center; L_seg_sin + y_center];
    Pos_X_2 = [-L_seg_cos + x_center; -L_seg_sin + y_center];
    Data_kept.Centroid_x(ind_mod(i)) = Pos_X_1(1,1);
    Data_kept.Centroid_y(ind_mod(i)) = Pos_X_1(2,1);
    Data_kept.axis_major_length(ind_mod(i)) = 2*L_seg;
    Line_to_add = Data_kept(ind_mod(i), :);
    Line_to_add.Mask_nb = New_ID(i); %New ID
    Line_to_add.Centroid_x = Pos_X_2(1,1);
    Line_to_add.Centroid_y = Pos_X_2(2,1);
    Data_kept = [Data_kept; Line_to_add]; % Addition new lines (order will be modified at the end of the main for loop)
 end
[~, sorted_time] = sort(Data_kept.Timepoint);
Data_kept = Data_kept(sorted_time, :);
end