function cost = weight_cost(Data_kept, Table_fin, weights)
    time_point = unique(Data_kept.Timepoint);
    cost = 0;
    for j = 1:length(time_point)
        time_max = time_point(j);
        ind_time_max = find(Data_kept.Timepoint == time_max);
        ID_stat_Data_kept = Data_kept.Mask_nb(ind_time_max);
        ind_time_max_2 = find(Table_fin.Timepoint == time_max);
        ID_stat_Data_kept_2 = Table_fin.Mother_ID(ind_time_max_2);
        % Lineage_1 = Data_kept.TrueLineage(ind_time_max); Lineage_2 = Table_fin.Lineage(ind_time_max_2);
        % cost = Lineage_2 - Lineage_1; cost(cost ~= 0) = 1;
        % cost = sum(cost);
        for i = 1:length(ID_stat_Data_kept)
            ind_ID = find(Table_fin.Mother_ID == ID_stat_Data_kept(i));
    %         Daughter_true = Data_kept.TrueLineage(ind_time_max(i)); 
            Daughter_true = [Data_kept.TrueID_1(ind_time_max(i)) Data_kept.TrueID_2(ind_time_max(i))];
    %         Daughter_pred = Table_fin.Lineage(ind_ID);
            Daughter_pred = [Table_fin.Daughter_ID_1(ind_ID) Table_fin.Daughter_ID_2(ind_ID)];
            matches = sum(ismember(Daughter_true, Daughter_pred));
            if matches == 2
                cost = cost + 0; % Both values are identical (order doesn't matter)
            elseif matches == 1
                cost = cost + 1; % One value is identical
            else
                cost = cost + 2; % Both values are different
            end
    %         if Daughter_true == Daughter_pred
    %             cost = cost + 0;
    %         else 
    %             cost = cost + 1;
    %         end
        end
    end
% regularization = 0.1*(sum(weights) - 1)^2; % Encourage weights to sum to 1
% cost = cost + regularization;
% cost = (weights(1) - 0.5)^2 + (weights(2) - 0.3)^2 + (weights(3) - 0.2)^2;
end