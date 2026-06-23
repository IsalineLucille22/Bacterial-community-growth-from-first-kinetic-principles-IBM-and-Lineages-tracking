function [Data_kept, Daughter_temp, ind_empty] = Disappearing_Cells(Data_kept, Daughter_temp, ind_empty, ind_t_1, t_2)
nb_empty = length(ind_empty);
if  ~isempty(ind_empty)
    New_ID = (max(Data_kept.Mask_nb) + 1):(max(Data_kept.Mask_nb) + nb_empty);
    ind_to_add = ind_t_1(ind_empty);%min(ind_t_1) + sum_h_13_row - 1;%
    Line_to_add = Data_kept(ind_to_add, :);
    Line_to_add.Timepoint = t_2*ones(nb_empty, 1);
    Line_to_add.Mask_nb = New_ID'; %New ID
    Line_to_add.KeptOrNot = Data_kept.KeptOrNot(ind_to_add) + 1;% zeros(nb_empty, 1);
    %%%%
    ind_temp = find(Line_to_add.KeptOrNot > 20); %Remove cells that are not present for more than x time points. If x = inf, never removed them.
    Line_to_add(ind_temp, :) = [];
    ind_empty(ind_temp) = [];
    nb_empty = nb_empty - sum(length(ind_temp));
    New_ID(ind_temp) = [];
    %%%%
    Data_kept = [Data_kept; Line_to_add]; % Addition new lines (order will be modified at the end of the main for loop)
    Daughter_temp(ind_empty, :) = New_ID'.*ones(nb_empty, 2); %Attribute their own copies as daugthers to make the bond with time t + 2
end
end