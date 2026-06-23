function Lineage_Mat_fin = fun_Lineage_Mat(Lineage_Mat)
Lineage_Mat_fin = [];
% nb_iter = length(Lineage_Mat);
max_ID = 0;
i = 1;
% while i < (nb_iter - 1) 
while ~isempty(Lineage_Mat{i + 1})
    temp_1 = Lineage_Mat{i};
    ID_1 = (max_ID + 1):(max_ID + length(temp_1));
    max_ID  = max(ID_1);
    lineage_temp = zeros(length(temp_1), 3);
    temp_2 = Lineage_Mat{i + 1};
    ID_2 = (max_ID + 1):(max_ID + length(temp_2));
    temp_2_low = temp_2(1:length(temp_1));
    temp_2_up = temp_2(length(temp_1) + 1:end);
    b = [];
    no_div = ismember(temp_1, temp_2_low); div = ~ismember(temp_1, temp_2_low);
    if ~isempty(temp_2_up)
        [~, b] = ismember(temp_2(div), temp_2_up - 1);
    end
    lineage_temp(no_div, 2) = ID_2(no_div);
    add_cell = ID_2((length(temp_1) + 1):end);
    lineage_temp(div, 2:3) = [ID_2(div); add_cell(b)]';
    lineage_temp(:, 1) = ID_1;
    Lineage_Mat_fin = [Lineage_Mat_fin; lineage_temp];
    i = i + 1;
end
% end
end