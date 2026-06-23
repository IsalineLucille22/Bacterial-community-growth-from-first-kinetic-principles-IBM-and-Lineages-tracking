function [Lineage_ind, val, Daughter_temp] = Link_New_Lineages(Dist, Lineage_ind, val, Daughter_temp, cell_ID_2, Threshold_val)
%%%% If some cells are not linked to others, check the distance to cells at time t−1 to see if there is a free "mother" cell that can be linked to them.
ind_new_Lineage = find(val == 0); %Find the new lineages
ind_mother_with_2_daugters = Daughter_temp(:, 1) - Daughter_temp(:, 2) ~= 0; %Find mothers with already 2 daughters (can't be used)
Dist_temp = Dist;
Dist_temp(ind_mother_with_2_daugters, :) = inf;
[val_min, ind_link] = min(Dist_temp(:, ind_new_Lineage), [], 1); %ind_link = index mother
% if length(ind_link) > 1
%     c = 10;
% end
unique_val = unique(ind_link);
for i = 1:length(unique_val)
    ind_link_temp = unique_val(i); %Potential mother
    ind_temp = find(ind_link == ind_link_temp);
    temp_val = val_min(ind_temp);
    potential_daughters = ind_new_Lineage(ind_temp);
    if Daughter_temp(ind_link_temp, 1)*Daughter_temp(ind_link_temp, 2) == 0
        [val_temp, temp] = mink(temp_val, 2); %The two minimal value
        ind_min_threshold = val_temp <= Threshold_val;
        val_temp = val_temp(ind_min_threshold); temp = temp(ind_min_threshold);
        if ~isempty(val_temp)
            new_daughters = potential_daughters(temp);
            Daughter_temp(ind_link_temp, :) = cell_ID_2(new_daughters);
            Lineage_ind(new_daughters) = ind_link_temp;
            val(new_daughters) = inf; %Test it
        end
    else
        [val_temp, temp] = min(temp_val); %The one minimal value
        if val_temp < 100
            new_daughters = potential_daughters(temp);
            Daughter_temp(ind_link_temp, 2) = cell_ID_2(new_daughters);
            Lineage_ind(new_daughters) = ind_link_temp;
            val(new_daughters) = inf; %Test it
        end
    end

% if val_min < inf
%     Lineage_ind(ind_new_Lineage) = ind_link;
%     if Daughter_temp(ind_link, 1).*Daughter_temp(ind_link, 2) == 0
%         Daughter_temp(ind_link, :) = cell_ID_2(ind_new_Lineage);
%     else
%         Daughter_temp(ind_link, 2) = cell_ID_2(ind_new_Lineage);
%     end
%     val(ind_new_Lineage) = inf; %Test it
% end
end
end