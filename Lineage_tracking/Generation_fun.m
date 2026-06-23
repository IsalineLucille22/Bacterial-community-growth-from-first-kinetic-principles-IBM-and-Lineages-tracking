function [Daughter_temp, generation_time, cell_ID_death, ind_Daughters] = Generation_fun(Daughter_temp, generation_time, cell_ID_death, cell_ID_1)
    generation_add = zeros(1, length(generation_time));
    ind_2_daught = Daughter_temp(:, 1) - Daughter_temp(:, 2) ~= 0; %Find indices with divison
    generation_add(ind_2_daught) = 1; %When division add + 1 to the generation time
    generation_time = generation_time + generation_add; %Generation time mother
    cell_ID_death(ind_2_daught) = cell_ID_1(ind_2_daught);
    ind_Daughters = [Daughter_temp(ind_2_daught, 1)' Daughter_temp(ind_2_daught, 2)']; %Keep the daughters indices
end
