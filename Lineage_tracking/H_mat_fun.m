function h_mat = H_mat_fun(Dist, n_1, n_2, vect_Cell_length_1, vect_Cell_length_2, height_cell_1, height_cell_2, Seg_tot_1, Seg_tot_2, weights, Seg_tot_1_rect, Seg_tot_2_rect, method_number)
Threshold_len = 2*(vect_Cell_length_1 + vect_Cell_length_2');
temp_kept = Dist <= Threshold_len;
h_mat = zeros(n_1, n_2);
% h_mat_2 = zeros(n_1, n_2);
for j = 1:n_1
    ind_temp = find(temp_kept(j, :) == 1); %Indices to test
    for k = 1:length(ind_temp)
        d_0_vect = (height_cell_1(j) + height_cell_2(ind_temp(k)))/2;%max(height_cell_1(j), height_cell_2(ind_temp(k)));%
        Over = zeros(1, 2);
        n_init = zeros(2, 2);
        X_o_init = {n_init, n_init};
        % h_mat(j, ind_temp(k)) = OverlapVal({Seg_tot_1{j}, Seg_tot_2{ind_temp(k)}}, d_0_vect, Over, n_init, X_o_init);
        dist_ind = Dist(j, ind_temp(k));
        h_mat(j, ind_temp(k)) = OverlapVal_v2({Seg_tot_1{j}, Seg_tot_2{ind_temp(k)}}, d_0_vect, Over, n_init, X_o_init, dist_ind, weights, height_cell_1(j), height_cell_2(ind_temp(k)), {Seg_tot_1_rect{j}, Seg_tot_2_rect{ind_temp(k)}}, method_number);
    end
end
end