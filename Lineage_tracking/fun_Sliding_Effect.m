function [Pos_X_1, Pos_X_2, cell_angle_1,  cell_angle_2, Dist_mat, Data_kept] = fun_Sliding_Effect(Data_kept, Pos_X_1, Pos_X_2, Dist_mat, cell_angle_1,  cell_angle_2, ind_t_1, ind_t_2)
[min_val, Lineage_ind_row] = min(Dist_mat, [], 1);
abs_min = mode(min_val); %min(min_val);
Lineage_ind_col = find(abs(min_val - abs_min) < 1e-08); Lineage_ind_col = Lineage_ind_col(1);
Lineage_ind_row = Lineage_ind_row(Lineage_ind_col);
if abs_min > 0 && abs_min < 20 %length(Lineage_ind_col) == 1 %&& abs_min < 10
    dist_x = Pos_X_1(1, Lineage_ind_row) - Pos_X_2(1, Lineage_ind_col); %x sliding
    dist_y = Pos_X_1(2, Lineage_ind_row) - Pos_X_2(2, Lineage_ind_col); %y sliding
    Pos_X_2(1, :) = Pos_X_2(1, :) + dist_x;
    Pos_X_2(2, :)  = Pos_X_2(2, :) + dist_y;
    Data_kept.Centroid_x(ind_t_2) = Pos_X_2(1, :);
    Data_kept.Centroid_y(ind_t_2) = Pos_X_2(2, :);
    % diff_orientation = cell_angle_1(Lineage_ind_row) - cell_angle_2(Lineage_ind_col);
    % cell_angle_1 = cell_angle_1 - diff_orientation;
    Dist_mat = distEuclid(Pos_X_1, Pos_X_2);
end 
% if abs_min < 60
    % unique_vals = unique(Lineage_ind_row);
    % [counts, edges] = histcounts(Lineage_ind_row, [unique_vals, max(unique_vals) + 1]); %Only if there is only positive values 
    % edges = edges(1:end - 1);
    % ind_edges = edges(counts == 1);
    % dist_x = 0; dist_y = 0; 
    % diff_orientation = 0;
    % for w = 1:length(ind_edges)
    %     ind_col = find(Lineage_ind_row == ind_edges(w));
    %     dist_x = dist_x + Pos_X_1(1, ind_edges(w)) - Pos_X_2(1, ind_col); %x sliding
    %     dist_y = dist_y + Pos_X_1(2, ind_edges(w)) - Pos_X_2(2, ind_col); %y sliding
    %     diff_orientation = diff_orientation + cell_angle_1(ind_edges(w)) - cell_angle_2(ind_col);
    % end
    % dist_x = dist_x/length(ind_edges);
    % dist_y = dist_y/length(ind_edges);
    % diff_orientation = diff_orientation/length(ind_edges);
    % Pos_X_1(1, :) = Pos_X_1(1, :) - dist_x;
    % Pos_X_1(2, :)  = Pos_X_1(2, :) - dist_y;
    % cell_angle_1 = cell_angle_1 - diff_orientation;
% end
end