function [Data_kept, Pos_X] = Fill_NaN_Pos(ind_t, Data_kept, Pos_X)
is_nan_ind = find(isnan(Pos_X(1,:)));
if sum(is_nan_ind) > 0
    min_ind = min(ind_t); ind_t = is_nan_ind  + min_ind - 1;
    Data_kept.Centroid_x(ind_t) = Data_kept.GeomX(ind_t);
    Pos_X(1,is_nan_ind) = Data_kept.GeomX(ind_t);
    Data_kept.Centroid_y(ind_t) = Data_kept.GeomY(ind_t);
    Pos_X(2,is_nan_ind) = Data_kept.GeomY(ind_t);
end
end