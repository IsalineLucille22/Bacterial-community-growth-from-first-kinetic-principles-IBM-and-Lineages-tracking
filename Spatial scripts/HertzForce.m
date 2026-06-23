function [h_mat, X_io_mat, X_jo_mat, n_ij_mat, X_io_grav] = HertzForce(Pos_X, i, vect_Cell_Length, ~, height_cell, Seg_tot, sp_nb, n, h_mat_init, X_io_mat_init, Over, n_init, X_o_init, X_io_grav)
r_1 = height_cell(sp_nb(i))/2;
% n = length(Pos_X(1,:));%Number of cells
h_mat = h_mat_init;%zeros(n, 1); %Put as agrugments?
% h_mat_2 = h_mat_init;
X_io_mat = X_io_mat_init;%zeros(n, 2);
X_jo_mat = X_io_mat;
n_ij_mat = X_io_mat;
r_2_vect = height_cell(sp_nb((i+1):n))/2;
d_0_vect = r_1 + r_2_vect;
dist = vecnorm([Pos_X(1,i) - Pos_X(1,(i+1):n); Pos_X(2,i) - Pos_X(2,(i+1):n)]); %Matlab version, less efficient
% dist = Dist_Mex(Pos_X(1,i) - Pos_X(1,(i+1):n), Pos_X(2,i) - Pos_X(2,(i+1):n));
threshold_forces = (vect_Cell_Length((i+1):n) + vect_Cell_Length(i))/2 + d_0_vect';
ind_2 = find(dist <= threshold_forces);
ind = ind_2 + i;
L = length(ind);
for k = 1:L
    % [h_mat(ind(k)), X_io_mat(ind(k), :), X_jo_mat(ind(k), :), n_ij_mat(ind(k), :)] = OverlapVal({Seg_tot{i}, Seg_tot{ind(k)}}, d_0_vect(ind_2(k)), Over, n_init, X_o_init);
    [h_mat(ind(k)), X_io_mat(ind(k), :), X_jo_mat(ind(k), :), n_ij_mat(ind(k), :)] = OverlapVal_v2({Seg_tot{i}, Seg_tot{ind(k)}}, d_0_vect(ind_2(k)), Over, n_init, X_o_init);
    % [h_mat_2(ind(k)), a, b, c] = OverlapVal_v2({Seg_tot{i}, Seg_tot{ind(k)}}, d_0_vect(ind_2(k)), Over, n_init, X_o_init);
    % if h_mat_2(ind(k)) ~= h_mat(ind(k))
    %     c = 10;
    % end
    %Only if gravitational force, be carful to consider ind(k) as the
    % seg_temp = Seg_tot{i};
    %bottom boundary.
    % if ind(k) == n - 1 %&& min(seg_temp(2,:)) <= (height_cell(sp_nb(i))/2)
    %     ind_temp = seg_temp(2,:) == min(seg_temp(2,:));
    %     X_io_grav = seg_temp(:,ind_temp);
    % end
end
% h_mat_temp = zeros(L, 1); X_io_mat_temp = zeros(L, 2); X_jo_mat_temp = zeros(L, 2);n_ij_mat_temp = zeros(L, 2);
% parfor k = 1:L
% %     [h_mat(ind(k)), X_io_mat(ind(k), :), X_jo_mat(ind(k), :), n_ij_mat(ind(k), :)] = OverlapVal({Seg_tot{i}, Seg_tot{ind(k)}}, d_0_vect(ind_2(k)));
%     [h_mat_temp(k), X_io_mat_temp(k, :), X_jo_mat_temp(k, :), n_ij_mat_temp(k, :)] = OverlapVal({Seg_tot{i}, Seg_tot{ind(k)}}, d_0_vect(ind_2(k)));
% %     [h_mat_temp, X_io_mat_temp, X_jo_mat_temp, n_ij_mat_temp] = OverlapVal({Seg_tot{i}, Seg_tot{ind}}, d_0_vect(ind_2));
% end
% if ~isempty(ind)
%     h_mat(ind) = h_mat_temp; X_io_mat(ind, :) = X_io_mat_temp;
%     X_jo_mat(ind, :) = X_jo_mat_temp; n_ij_mat(ind, :) = n_ij_mat_temp;
% end
