function [h_mat, Lineage_1, Lineage_2, Lineage_ind, Daughter_temp, nb_Cell_t_1, nb_Cell_t_2, Data_kept, cell_ID_1, cell_ID_2, Kept_NOR_1, Kept_NOR_2, height_cell_1, height_cell_2, vect_Cell_length_1, vect_Cell_length_2, cell_angle_1, cell_angle_2, Pos_X_1, Pos_X_2] = Replace_big_Cells(Lineage_1, i, time_point, Data_kept, Lineage_ind, h_mat, ind_empty, height_cell_1, height_cell_2, vect_Cell_length_1, vect_Cell_length_2, cell_ID_2, ind_t_1, ind_t_2, t_2, real_sim, overlap_val_margin, max_Lineage)
sum_row = sum(h_mat, 2);
moth_no_daughter = intersect(find(sum_row > 0), ind_empty);
% ind_Mother_used = [];
for s = 1:length(moth_no_daughter)
    New_ID = (max(Data_kept.Mask_nb) + 1);
    [~, ind_daughter] = max(h_mat(moth_no_daughter(s), :));
    ind_true_mother = Lineage_ind(ind_daughter);
    if (height_cell_2(ind_daughter) > 1.5*height_cell_1(ind_true_mother) || vect_Cell_length_2(ind_daughter) > 1.5*vect_Cell_length_1(ind_true_mother))  && (h_mat(moth_no_daughter(s), ind_daughter) > 0.5*h_mat(ind_true_mother, ind_daughter))
        ID_daugther = cell_ID_2(ind_daughter);
        ind_daughter = ind_t_2(ind_daughter);
        % is_new_mother = intersect(ind_true_mother, ind_Mother_used);
        % if  ~isempty(is_new_mother)
        %     c = 10;
        % end
        ind_true_mother = ind_t_1(ind_true_mother);
        Data_kept(ind_daughter, :) = Data_kept(ind_true_mother, :);
        Data_kept.Mask_nb(ind_daughter) = ID_daugther;
        Data_kept.Timepoint(ind_daughter) = t_2;
        ind_to_add = ind_t_1(moth_no_daughter(s));
        Line_to_add = Data_kept(ind_to_add, :);
        Line_to_add.Timepoint = t_2;
        Line_to_add.Mask_nb = New_ID'; %New ID
        Line_to_add.KeptOrNot = 1;
        Data_kept = [Data_kept; Line_to_add]; % Addition new lines (order will be modified at the end of the main for loop)
        % ind_Mother_used = [ind_Mother_used ind_true_mother]; 
    end
end

[~, sorted_time] = sort(Data_kept.Timepoint);
Data_kept = Data_kept(sorted_time, :);

t_1 = time_point(i);
t_2 = time_point(i + 1);

ind_t_1 = find(Data_kept.Timepoint == t_1);
nb_Cell_t_1 = length(ind_t_1);
Data_1 = Data_kept(ind_t_1, :);
ind_t_2 = find(Data_kept.Timepoint == t_2);
nb_Cell_t_2 = length(ind_t_2);
Data_2 = Data_kept(ind_t_2, :);
Pos_X_1 = [Data_1.Centroid_x, Data_1.Centroid_y]';
Pos_X_2 = [Data_2.Centroid_x, Data_2.Centroid_y]';
Dist = distEuclid(Pos_X_1, Pos_X_2);

vect_Cell_length_1 = Data_1.axis_major_length; 
height_cell_1 = Data_1.axis_minor_length;
cell_angle_1 = Data_1.orientation;
Kept_NOR_1 = Data_1.KeptOrNot;
vect_Cell_length_2 = Data_2.axis_major_length; 
height_cell_2 = Data_2.axis_minor_length;
cell_angle_2 = Data_2.orientation;
Kept_NOR_2 = Data_2.KeptOrNot;

vect_Cell_length_1 = max(vect_Cell_length_1, mean(Data_kept.axis_major_length(ind_t_1)));
n_1 = height(Data_1(:, 1));
Seg_tot_1 = arrayfun(@(x) Rect2Seg([Pos_X_1(1,x) Pos_X_1(2,x) (vect_Cell_length_1(x) - real_sim*height_cell_1(x)) height_cell_1(x)], cell_angle_1(x)),1:n_1,'UniformOutput',false); %Find an alternative
cell_ID_1 = Data_1.Mask_nb;
height_cell_1 = Data_kept.axis_minor_length(ind_t_1);%mean(Data_kept.axis_minor_length(ind_t_1))*ones(nb_Cell_t_1, 1);

vect_Cell_length_2 = max(vect_Cell_length_2, mean(Data_kept.axis_major_length(ind_t_2)));
n_2 = height(Data_2(:, 1)); 
Seg_tot_2 = arrayfun(@(x) Rect2Seg([Pos_X_2(1,x) Pos_X_2(2,x) (vect_Cell_length_2(x) - real_sim*height_cell_2(x)) height_cell_2(x)], cell_angle_2(x)),1:n_2,'UniformOutput',false); %Find an alternative
cell_ID_2 = Data_2.Mask_nb;
height_cell_2 = Data_kept.axis_minor_length(ind_t_2);%mean(Data_kept.axis_minor_length(ind_t_2))*ones(nb_Cell_t_2, 1);

Daughter_temp = zeros(n_1, 2);

h_mat = H_mat_fun(Dist, n_1, n_2, vect_Cell_length_1 - real_sim*height_cell_1, vect_Cell_length_2 - real_sim*height_cell_2, overlap_val_margin*height_cell_1, overlap_val_margin*height_cell_2, Seg_tot_1, Seg_tot_2); %H_mat_fun(Dist, n_1, n_2, vect_Cell_length_1 - height_cell_1, vect_Cell_length_2 - height_cell_2, overlap_val_margin*height_cell_1, overlap_val_margin*height_cell_2, Seg_tot_1, Seg_tot_2);
h_mat(h_mat < 0) = 0;

[val, Lineage_ind] = max(h_mat, [], 1);
Lineage_2 = Lineage_1(Lineage_ind);
nb_zeros = sum(val == 0);
Lineage_2(val == 0) = max_Lineage + 1:max_Lineage + nb_zeros;
Lineage_ind(val == 0) = 0;
end