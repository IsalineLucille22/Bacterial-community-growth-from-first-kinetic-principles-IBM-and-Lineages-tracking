function [Data_kept, vect_Cell_length, height_cell, cell_angle] = Fill_NaN(ind_t, Dist, Data_kept, vect_Cell_length, height_cell, cell_angle, Pos_X, vect_Cell_length_opp, height_cell_opp, cell_angle_opp, area_opp)
is_nan_ind = find(isnan(cell_angle));
if sum(is_nan_ind) > 0
    is_nan_ind_opp = isnan(cell_angle_opp);
    Dist(is_nan_ind_opp, :) = inf; 
    ind_t = ind_t(is_nan_ind);
    [~, ind_min] = min(Dist, [], 1);
    vect_Cell_length(is_nan_ind) = vect_Cell_length_opp(ind_min(is_nan_ind));
    Data_kept.Centroid_x(ind_t) = Pos_X(1, is_nan_ind); %Data_kept.centroid_1(ind_t) = Pos_X(1, is_nan_ind);
    Data_kept.Centroid_y(ind_t) = Pos_X(2, is_nan_ind); %Data_kept.centroid_0(ind_t) = Pos_X(2, is_nan_ind);
    Data_kept.axis_major_length(ind_t) = vect_Cell_length_opp(ind_min(is_nan_ind));
    Data_kept.axis_minor_length(ind_t) = height_cell_opp(ind_min(is_nan_ind));
    Data_kept.orientation(ind_t) = cell_angle_opp(ind_min(is_nan_ind));%Data_kept.orientation_radian_(ind_t) = cell_angle_opp(ind_min(is_nan_ind));
    if ~isempty(area_opp)
        Data_kept.area(ind_t) = area_opp(ind_min(is_nan_ind));
    end
    height_cell(is_nan_ind) = height_cell_opp(ind_min(is_nan_ind));
    cell_angle(is_nan_ind) = cell_angle_opp(ind_min(is_nan_ind));
end
end