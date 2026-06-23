function [h_mat, Daughter_temp, Life_time, Life_add, ind_empty] = Life_Time_fun_v2(h_mat, Daughter_temp, Life_time, Lineage_ind, cell_ID_2)
Life_add = ones(1, length(Life_time));
ind_2_daught = Daughter_temp(:, 1) - Daughter_temp(:, 2) ~= 0;
Life_add(ind_2_daught) = 0;
Life_time_temp = Life_time + Life_add; %Life time mother
mother_twice_div = find(Life_time_temp == 0); 
col_daught = zeros(length(Life_time), length(mother_twice_div));
ind_daughter_to_rem_tot = []; ind_daughter_to_kept_tot = [];
%Reallocate the cells that divide twice to another branch
for i = 1:length(mother_twice_div) %min(1, length(mother_twice_div))
    h_mat_temp = zeros(size(h_mat));
    ind_daughter = find(Lineage_ind == mother_twice_div(i));
    %%%%
    linear_indices = sub2ind(size(h_mat), Lineage_ind(ind_daughter), ind_daughter); 
    [~, min_ind] = min(h_mat(linear_indices));
    %%%%
    h_mat_temp((Life_time_temp > 1) & (Life_add > 0), :) = h_mat((Life_time_temp > 1) & (Life_add > 0), :); %Cell whose last divison > 1 and that doesn't divide 
    ind_daughter_to_rem = ind_daughter(min_ind);
    col_daught(:, i) = h_mat_temp(:, ind_daughter_to_rem); %Here only keep the mothers from the same lineage
    ind_daughter_to_rem_tot = [ind_daughter_to_rem_tot, ind_daughter_to_rem];
    ind_daughter_to_kept_tot = [ind_daughter_to_kept_tot, ind_daughter(1:2 ~= min_ind)];
end
id_b = 1:length(Life_time);
[~, id_a] = max(col_daught, [], 2);
if ~isempty(id_a)
    linear_indices = sub2ind(size(col_daught), id_b', id_a); 
    col_daught_old = col_daught;
    col_daught = zeros(length(Life_time), length(mother_twice_div));
    col_daught(linear_indices) = col_daught_old(linear_indices);
    h_mat_temp(:, ind_daughter_to_rem_tot) = col_daught;
end

%ind_daugther, min_ind
for i = 1:length(mother_twice_div) 
    ind_daughter_to_rem = ind_daughter_to_rem_tot(i);
    [val, pot_mother_ind] = max(col_daught(:, i));
    if val > 0
        if Daughter_temp(pot_mother_ind, 1) == 0
            Daughter_temp(pot_mother_ind, :) = cell_ID_2(ind_daughter_to_rem)*ones(1,2); %If zeros no daughter so cell_ID_2(ind_daughter_to_rem) because the unique daughter
        else
            Daughter_temp(pot_mother_ind, 2) = cell_ID_2(ind_daughter_to_rem); %cell_ID_2(ind_daughter_to_rem) because the second daughter 
        end
        Daughter_temp(mother_twice_div(i), :) = cell_ID_2(ind_daughter_to_kept_tot(i)); %Only keep the ID of the remaining daughter
        h_mat(:, ind_daughter_to_rem) = h_mat_temp(:, ind_daughter_to_rem);
    end
end
Life_add = ones(1, length(Life_time));
ind_empty = find(Daughter_temp(:,1).*Daughter_temp(:,2) == 0);
ind_2_daught = Daughter_temp(:, 1) - Daughter_temp(:, 2) ~= 0;
Life_add(ind_2_daught) = 0;
Life_time = Life_time + Life_add; %Life time mother
end