close all
plot([A(1), B(1)], [A(2), B(2)], 'b--', 'LineWidth', 2); 
hold on;
plot([pt_1(1) pt_2(1)], [pt_1(2) pt_2(2)], 'r-.', 'LineWidth', 1);
hold on;
plot([proj_1(1) pt_1(1)], [proj_1(2) pt_1(2)], '--x', 'LineWidth', 0.5);
hold on;
plot([proj_2(1) pt_2(1)], [proj_2(2) pt_2(2)], '--x', 'LineWidth', 0.5);
hold on;
plot(pt_inter(1), pt_inter(2), 'o', 'LineWidth', 6);

close all
plot([A(1), B(1)], [A(2), B(2)], 'b--', 'LineWidth', 2); 
hold on;
plot([pt_1(1) pt_2(1)], [pt_1(2) pt_2(2)], 'r-.', 'LineWidth', 1);
hold on;
plot([X_io(1) X_jo(1)], [X_io(2) X_jo(2)], '--k', 'LineWidth', 2);
axis equal

close all
plot([A(1), B(1)], [A(2), B(2)], 'b--', 'LineWidth', 2); 
hold on;
plot([pt_1(1) pt_2(1)], [pt_1(2) pt_2(2)], 'r-.', 'LineWidth', 1);
hold on;
plot([proj_C_on_AB(1) pt_1(1)], [proj_C_on_AB(2) pt_1(2)], '--x', 'LineWidth', 0.5);
hold on;
plot([proj_D_on_AB(1) pt_2(1)], [proj_D_on_AB(2) pt_2(2)], '--x', 'LineWidth', 0.5);
axis equal

unique_vals = unique(Lineage_ind(Lineage_ind ~= 0));
h_mat_old = h_mat;
if ~isempty(unique_vals)
    counts = histcounts(Lineage_ind, [unique_vals, max(unique_vals) + 1]);
    for j = 1:n_1 
        ind_temp = find(Lineage_ind == j);
        ID_Daughters = cell_ID_2(ind_temp);
        if length(ID_Daughters) > 2
            [~, b] = maxk(h_mat(j, ind_temp), 2);
            Ind_Daughters_NK = ind_temp(setdiff(1:length(ID_Daughters), b));
            h_mat(:, Ind_Daughters_NK) = 0;
        end
    end
    [val_2_daught, Lineage_ind_2_daught] = max(h_mat, [], 1);
    Lineage_ind_2_daught(val_2_daught == 0) = 0;
    ind_daughter = 1:s;
    ind_daughter(val_2_daught == 0) = [];
    linear_indices_2_daught = sub2ind(size(h_mat_old), Lineage_ind_2_daught(Lineage_ind_2_daught ~= 0), ind_daughter);
end

h_mat(linear_indices_2_daught) = inf;