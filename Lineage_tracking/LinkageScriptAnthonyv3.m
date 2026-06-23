clear
close all
clc

%Save or Not
save_data = 1; %1 if the video is saved, 0 otherwise
%The final tables, cellcounts and Lineage variables are automatically saved
%in the folder Data. If their names are not changed, it will overwrite the
%existing variables.

%Loadind data
%Anthony's data
% Data_kept_init = readtable(strcat('Data/','final_merged_table_pos1.xlsx'), 'Sheet', 1,'Format','auto'); 
Data_kept_init = readtable(strcat('Data/','final_merged_table_Pos_double.xlsx'), 'Sheet', 1,'Format','auto'); 
% Data_kept_init = readtable('/Volumes/RECHERCHE/FAC/FBM/DMF/jvanderm/ice-clc/D2c/RawData/Anthony_C/Microscopy/Timelapse/230124_7639_ahpc/_1/Trimmed/_1/Pos11/final_merged_table.xlsx');
% Data_kept_init = readtable('/Volumes/RECHERCHE/FAC/FBM/DMF/jvanderm/ice-clc/D2c/RawData/Anthony_C/Microscopy/Timelapse/221124_7639/_1/Pos16/final_merged_table.xlsx');
% Data_kept_init = readtable('/Volumes/RECHERCHE/FAC/FBM/DMF/jvanderm/ice-clc/D2c/RawData/Anthony_C/Microscopy/Timelapse/230124_7639_ahpc/_1/Trimmed/_1/Pos2/final_merged_table.xlsx');
% Data_kept_init = readtable('/Volumes/RECHERCHE/FAC/FBM/DMF/jvanderm/ice-clc/D2c/RawData/Anthony_C/Microscopy/Timelapse/231206_07fpr/_1/Pos1/final_merged_table.xlsx'); %Heart-shaped
% Data_kept_init = readtable('/Volumes/RECHERCHE/FAC/FBM/DMF/jvanderm/ice-clc/D2c/RawData/Anthony_C/Microscopy/Timelapse/231026_1683_3317/_1/Pos4/final_merged_table.xlsx');
% Data_kept_init = readtable('/Volumes/RECHERCHE/FAC/FBM/DMF/jvanderm/ice-clc/D2c/RawData/Anthony_C/Microscopy/Timelapse/231026_1683_3317/_1/Pos3/final_merged_table.xlsx');
num_sheet = 'Anthony_Pos1';%'merged_table_202312_Rahnella_Mono';%'Anthony_Pos1';%3;%'20230809_Curto_Mono_Pos8.xlsx';%'Anthony_Pos1';%'AnthonyPos1Corr';%2;%'BoukePos9Corr';%2; %'Pos14';%3;

Data_kept_init.Mask_nb = (1:height(Data_kept_init))'; %reinitialize cell numbers
Data_kept_init.Timepoint = Data_kept_init.Timepoint - min(Data_kept_init.Timepoint); %Modify time to start at 0 anyway

real_sim = 1; %1 if real, 0 if sim
Data_kept_init.Timepoint = Data_kept_init.Timepoint - min(Data_kept_init.Timepoint); %Modify time to start at 0 anyway
[~, sorted_time] = sort(Data_kept_init.Timepoint);
Data_kept_init = Data_kept_init(sorted_time, :);
Data_kept_init.orientation = -(pi/2 - Data_kept_init.orientation);
Data_kept_init.Centroid_y = 2050 - Data_kept_init.Centroid_y;
Data_kept_init = Arrage_Timepoint_fun(Data_kept_init); %Sometimes, there are missing time points, which can cause issues later. To address this, we force time points to increment only by 1.
time_point = unique(Data_kept_init.Timepoint);
time_point_init = 1;
Stat_t = Stat_Phase_fun(Data_kept_init);
Data_kept_init = fun_Slinding_Effect_v2(Data_kept_init);
weights = [0.25, 0.25, 0.25, 0.25]; %Not used
method_number = 0; %If the new overlapping method is used, 1 otherwise
if method_number == 0
    real_sim = -1; %1 if real, 0 if sim
end

len_time_point = Stat_t;
overlap_val_margin_tot = 1:12; %Sequence for the overlapping margings (intentional increase of the cell diameter)
nb_No_unique = inf*ones(1, length(overlap_val_margin_tot)); Table_fin_tot = {};
int_num_cells = length(find(Data_kept_init.Timepoint == min(time_point)));
new_column = ones(height(Data_kept_init), 1);  %Generate new data for the column
Data_kept_init.KeptOrNot = new_column;  %Assign the new column
new_column = zeros(height(Data_kept_init), 1);  %Generate new data for the column
Data_kept_init.Life_time = new_column;  %Assign the new column

mean_length =  mean(Data_kept_init.axis_major_length(~isnan(Data_kept_init.axis_major_length))); %mean(Data_1.axis_major_length(~isnan(Data_1.axis_major_length)));%
mean_d =  mean(Data_kept_init.axis_minor_length(~isnan(Data_kept_init.axis_minor_length)));%mean(Data_1.axis_minor_length(~isnan(Data_1.axis_minor_length)));%
%%%% Comment or not this part if cell shape corrections should be applied 
Data_kept_init = Remove_Cells(Data_kept_init, 'axis_minor_length', mean_d, 0.2, 3); %Remove cells with a too small/big length
Data_kept_init = Length_cells_mod(Data_kept_init, mean_length, 4); %Split into 2 cells the cells that are too long 
Data_kept_init = Correction_Cell_Sizes(Data_kept_init, 'axis_minor_length', mean_d, 0.8, 1.5); %Correction of cell diameters. Modify the cell diameters that are too small/big in comparison to the average on the entire data set
Data_kept_init = Correction_Cell_Sizes(Data_kept_init, 'axis_major_length', mean_length, 0.4, 2.5); %Correction of cell lengths. Modify the cell lengths that are too small/big in comparison to the average on the entire data set
ind_t_1 = find(Data_kept_init.Timepoint == time_point(time_point_init));
nb_Cell_t_1 = length(ind_t_1);
Data_1 = Data_kept_init(ind_t_1, :);
Tab_fin  = [];
Data_kept = Data_kept_init;
n_1 = sum(Data_kept.Timepoint == time_point(time_point_init));
Lineage_2 = 1:n_1;
To_add_next = {};
max_Lineage = max(Lineage_2);

%%%%
Life_time = zeros(1, nb_Cell_t_1);
Life_time_daughter = Life_time;
generation_time = zeros(1, nb_Cell_t_1);
cell_ID_death = zeros(1, nb_Cell_t_1);
cell_ID_birth = Data_1.Mask_nb';
cell_ID_plot = ones(1, nb_Cell_t_1);
%%%%
Threshold_val = 100; %Maximum distance to be attributed to an existing lineage
for i = time_point_init:(len_time_point - 2)
    t_1 = time_point(i);
    t_2 = time_point(i + 1);
    t_3 = time_point(i + 2);

    ind_t_1 = find(Data_kept.Timepoint == t_1);
    nb_Cell_t_1 = length(ind_t_1);
    Data_1 = Data_kept(ind_t_1, :);
    ind_t_2 = find(Data_kept.Timepoint == t_2);
    nb_Cell_t_2 = length(ind_t_2);
    ind_t_3 = find(Data_kept.Timepoint == t_3);
    Data_3 = Data_kept(ind_t_3, :);
    Data_2 = Data_kept(ind_t_2, :);
    Pos_X_1 = [Data_1.Centroid_x, Data_1.Centroid_y]';
    [Data_kept, Pos_X_1] = Fill_NaN_Pos(ind_t_1, Data_kept, Pos_X_1);%Only for Bouke's data because of the different structure
    Pos_X_2 = [Data_2.Centroid_x, Data_2.Centroid_y]';
    [Data_kept, Pos_X_2] = Fill_NaN_Pos(ind_t_2, Data_kept, Pos_X_2);%Only for Bouke's data because of the different structure
    Pos_X_3 = [Data_3.Centroid_x, Data_3.Centroid_y]';
    [Data_kept, Pos_X_3] = Fill_NaN_Pos(ind_t_3, Data_kept, Pos_X_3);%Only for Bouke's data because of the different structure
    Dist = distEuclid(Pos_X_1, Pos_X_2);
    Dist_23 = distEuclid(Pos_X_2, Pos_X_3);
    Dist_23 = [Dist' Dist_23];

    vect_Cell_length_1 = Data_1.axis_major_length; 
    height_cell_1 = Data_1.axis_minor_length;
    cell_angle_1 = Data_1.orientation;
    Kept_NOR_1 = Data_1.KeptOrNot;
    GFP_cell_1 = Data_1.GFP_mean; %Also fill NaN values?
    MCHE_cell_1 = Data_1.mCherry_mean; %Also fill NaN values?
    vect_Cell_length_2 = Data_2.axis_major_length; 
    height_cell_2 = Data_2.axis_minor_length;
    cell_angle_2 = Data_2.orientation;
    Kept_NOR_2 = Data_2.KeptOrNot;
    vect_Cell_length_3 = Data_3.axis_major_length; 
    height_cell_3 = Data_3.axis_minor_length;
    cell_angle_3 = Data_3.orientation;
    GFP_cell_2 = Data_2.GFP_mean; %Also fill NaN values?
    MCHE_cell_2 = Data_2.mCherry_mean; %Also fill NaN values?

    %Fill the NaN values
    %Combination t-1 and t+1 for t_2
    [Data_kept, vect_Cell_length_1, height_cell_1, cell_angle_1] = Fill_NaN(ind_t_1, Dist', Data_kept, vect_Cell_length_1, height_cell_1, cell_angle_1, Pos_X_1, vect_Cell_length_2, height_cell_2, cell_angle_2, []);
    [Data_kept, vect_Cell_length_2, height_cell_2, cell_angle_2] = Fill_NaN(ind_t_2, Dist_23', Data_kept, vect_Cell_length_2, height_cell_2, cell_angle_2, Pos_X_2, [vect_Cell_length_1; vect_Cell_length_3], [height_cell_1; height_cell_3], [cell_angle_1; cell_angle_3], []);

    %Put here the sliding function and re-compute position updated %Comment
    %for Bouke
    % [Pos_X_1, Pos_X_2,  cell_angle_1,  cell_angle_2, Dist] = fun_Sliding_Effect(Pos_X_1, Pos_X_2, Dist, cell_angle_1,  cell_angle_2);

    vect_Cell_length_1 = max(vect_Cell_length_1, mean(Data_kept.axis_major_length(ind_t_1)));
    n_1 = nb_Cell_t_1; %height(Data_1(:, 1));
    Seg_tot_1 = arrayfun(@(x) Rect2Seg([Pos_X_1(1,x) Pos_X_1(2,x) (vect_Cell_length_1(x) - real_sim*height_cell_1(x)) height_cell_1(x)], cell_angle_1(x)),1:n_1,'UniformOutput',false); %Find an alternative
    Seg_tot_1_rect = arrayfun(@(x) Rect2Seg([Pos_X_1(1,x) Pos_X_1(2,x) (vect_Cell_length_1(x) - real_sim*height_cell_1(x)) height_cell_1(x)], cell_angle_1(x)),1:n_1,'UniformOutput',false); %Find an alternative
    cell_ID_1 = Data_1.Mask_nb;
    height_cell_1 = Data_kept.axis_minor_length(ind_t_1);%mean(Data_kept.axis_minor_length(ind_t_1))*ones(nb_Cell_t_1, 1);
    
    vect_Cell_length_2 = max(vect_Cell_length_2, mean(Data_kept.axis_major_length(ind_t_2)));
    n_2 = nb_Cell_t_2; %height(Data_2(:, 1)); 
    Seg_tot_2 = arrayfun(@(x) Rect2Seg([Pos_X_2(1,x) Pos_X_2(2,x) (vect_Cell_length_2(x) - real_sim*height_cell_2(x)) height_cell_2(x)], cell_angle_2(x)),1:n_2,'UniformOutput',false); %Find an alternative
    Seg_tot_2_rect = arrayfun(@(x) Rect2Seg([Pos_X_2(1,x) Pos_X_2(2,x) (vect_Cell_length_2(x) - real_sim*height_cell_2(x)) height_cell_2(x)], cell_angle_2(x)),1:n_2,'UniformOutput',false); %Find an alternative
    cell_ID_2 = Data_2.Mask_nb;
    height_cell_2 = Data_kept.axis_minor_length(ind_t_2);%mean(Data_kept.axis_minor_length(ind_t_2))*ones(nb_Cell_t_2, 1);
    
    m = 1; val_temp = 0;
    Lineage_1 = Lineage_2;
    Life_time = Life_time_daughter;
    Daughter_temp = zeros(n_1, 2);
    linear_indices = [];
    h_mat = H_mat_fun(Dist, n_1, n_2, vect_Cell_length_1 - real_sim*height_cell_1, vect_Cell_length_2 - real_sim*height_cell_2, height_cell_1, height_cell_2, Seg_tot_1, Seg_tot_2, weights, Seg_tot_1, Seg_tot_2, method_number ); 
    h_mat(h_mat < 0) = 0;
    h_mat(linear_indices) = inf;
    [val, Lineage_ind] = max(h_mat, [], 1);
    while m < length(overlap_val_margin_tot) && sum(val_temp == 0) > 0
    
        while isequal((val == 0), (val_temp == 0)) && m < length(overlap_val_margin_tot)
            overlap_val_margin = overlap_val_margin_tot(m);
            h_mat = H_mat_fun(Dist, n_1, n_2, vect_Cell_length_1 - real_sim*height_cell_1, vect_Cell_length_2 - real_sim*height_cell_2, overlap_val_margin*height_cell_1, overlap_val_margin*height_cell_2, Seg_tot_1, Seg_tot_2, weights, Seg_tot_1_rect, Seg_tot_2_rect, method_number ); 
            h_mat(h_mat < 0) = 0;
            h_mat(linear_indices) = inf;
            [val, Lineage_ind] = max(h_mat, [], 1);
            m = m + 1;
        end

        Lineage_2 = Lineage_1(Lineage_ind);
        nb_zeros = sum(val == 0);
        Lineage_2(val == 0) = max_Lineage + 1:max_Lineage + nb_zeros;
        Lineage_ind(val == 0) = 0;
        ind_empty = [];
    
        %%%% To add or remove. Impose unique lineage value for each tree
        unique_vals = unique(Lineage_ind(Lineage_ind ~= 0));
        if ~isempty(unique_vals)
            counts = histcounts(Lineage_ind, [unique_vals, max(unique_vals) + 1]);
            while max(counts) > 2
                for j = 1:n_1 
                    ind_temp = find(Lineage_ind == j);
                    ID_Daughters = cell_ID_2(ind_temp);
                    if length(ID_Daughters) > 2
                        [~, b] = maxk(h_mat(j, ind_temp), 2);
                        Ind_Daughters_NK = ind_temp(setdiff(1:length(ID_Daughters), b));
                        h_mat(j, Ind_Daughters_NK) = 0;
                        ID_Daughters = ID_Daughters(b);
                    end
                    if ~isempty(ID_Daughters)
                        Daughter_temp(j, :) = ID_Daughters;
                    end
                end
                [val, Lineage_ind] = max(h_mat, [], 1);
                Lineage_2 = Lineage_1(Lineage_ind);
                nb_zeros = sum(val == 0);
                Lineage_2(val == 0) = max_Lineage + 1:max_Lineage + nb_zeros;
                Lineage_ind(val == 0) = 0;
                unique_vals = unique(Lineage_ind(Lineage_ind ~= 0));
                counts = histcounts(Lineage_ind, [unique_vals, max(unique_vals) + 1]);
            end
            for j = 1:n_1 
                ind_temp = find(Lineage_ind == j);
                ID_Daughters = cell_ID_2(ind_temp);
                if length(ID_Daughters) > 2
                    [~, b] = maxk(h_mat(j, ind_temp), 2);
                    Ind_Daughters_NK = ind_temp(setdiff(1:length(ID_Daughters), b));
                    h_mat(j, Ind_Daughters_NK) = 0;
                    ID_Daughters = ID_Daughters(b);
                end
                if ~isempty(ID_Daughters)
                    Daughter_temp(j, :) = ID_Daughters;
                end
            end
            ind_empty = find(Daughter_temp(:,1).*Daughter_temp(:,2) == 0);%Remove
            [val, Lineage_ind] = max(h_mat, [], 1);
            Lineage_2 = Lineage_1(Lineage_ind);
            nb_zeros = sum(val == 0);
            Lineage_2(val == 0) = max_Lineage + 1:max_Lineage + nb_zeros;
            Lineage_ind(val == 0) = 0;
        end
        %%%%%
        ind_daughter = 1:n_2;
        ind_daughter(val == 0) = [];
        if ~isempty(ind_daughter)
            linear_indices = sub2ind(size(h_mat), Lineage_ind(Lineage_ind ~= 0), ind_daughter);  %Replace by inf the indices corresponding to a relation mother-daughter attributed
        end
        m = m + 1;
        val_temp = val;
    end

    %%%
    Daughter_temp_1 = Daughter_temp;
    Daughter_temp_1(Daughter_temp(:, 1) - Daughter_temp(:, 2) == 0, 2) = 0;
    Daughters_vect = [Daughter_temp_1(:,1); Daughter_temp_1(:,2)];
    Daughters_vect = Daughters_vect(Daughters_vect ~= 0);
    unique_vals_Daughters = unique(Daughters_vect);
    more_two_daughters = histcounts(Daughters_vect', [unique_vals_Daughters', max(unique_vals_Daughters) + 1]);
    if max(more_two_daughters) > 1
        c = 10;
    end
    %%%

    %%%% Life_time
    [h_mat, Daughter_temp, Life_time, Life_add, ind_empty] = Life_Time_fun_v2(h_mat, Daughter_temp, Life_time, Lineage_ind, cell_ID_2);
    Life_time_temp = Life_time;
    Life_time_temp(Life_add == 0) = 0;
    [val, Lineage_ind] = max(h_mat, [], 1);
    Life_time_daughter = Life_time_temp(Lineage_ind);
    Life_time_daughter(val == 0) = 0;
    Lineage_2 = Lineage_1(Lineage_ind);
    nb_zeros = sum(val == 0); 
    %%%%
    if nb_zeros > 0
        [Lineage_ind, val, Daughter_temp] = Link_New_Lineages(Dist, Lineage_ind, val, Daughter_temp, cell_ID_2, Threshold_val);
        ind_empty = find(Daughter_temp(:,1).*Daughter_temp(:,2) == 0);
        Lineage_2 = Lineage_1(Lineage_ind);
        nb_zeros = sum(val == 0); 
    end
    %%%%
    Lineage_2(val == 0) = max_Lineage + 1:max_Lineage + nb_zeros; %New lineages. Look at this and try to connect it to existing lineages by minimal distance
    %%%%Addition of disappearing cells
    [Data_kept, Daughter_temp, ind_empty] = Disappearing_Cells(Data_kept, Daughter_temp, ind_empty, ind_t_1, t_2);
    %%%%
    %%%%
    generation_time_save = generation_time;
    [Daughter_temp, generation_time, cell_ID_death, ind_Daughters] = Generation_fun(Daughter_temp, generation_time, cell_ID_death, cell_ID_1);
    generation_time_daughter = Daughters_Update_fun(generation_time, val, Lineage_ind);
    cell_ID_death_daughter = Daughters_Update_fun(cell_ID_death, val, Lineage_ind);
    cell_ID_birth_daughter = Daughters_Update_fun(cell_ID_birth, val, Lineage_ind);
    cell_ID_birth_daughter(val == 0) = cell_ID_2(val == 0);
    ind_div = ismember(cell_ID_2, ind_Daughters);
    cell_ID_birth_daughter(ind_div) = cell_ID_2(ind_div);
    cell_ID_plot_daughter = Daughters_Update_fun(cell_ID_plot, val, Lineage_ind);
    cell_ID_plot_daughter(val == 0) = 1;
    cell_ID_plot_daughter = plot_table_fun(cell_ID_plot_daughter, cell_ID_2, ind_Daughters, generation_time_daughter);
    %%%%
    Daughter_temp(Daughter_temp(:, 1) - Daughter_temp(:, 2) == 0, 2) = 0;
    Tab_1 = [cell_ID_1, Daughter_temp, Lineage_1', Data_1.Timepoint,  Pos_X_1', cell_angle_1, height_cell_1, vect_Cell_length_1, Kept_NOR_1, Life_time', generation_time_save', cell_ID_death', cell_ID_birth', GFP_cell_1, MCHE_cell_1, cell_ID_plot']; %Create here the all_tracks table 
    Tab_fin = [Tab_fin; Tab_1];
    max_Lineage = max_Lineage + nb_zeros;
    [~, sorted_time] = sort(Data_kept.Timepoint);
    Data_kept = Data_kept(sorted_time, :);
    Lineage_2 = [Lineage_2 Lineage_1(ind_empty)];
    Life_time_daughter = [Life_time_daughter Life_time(ind_empty)];
    %%%%%
    generation_time = [generation_time_daughter generation_time(ind_empty)];
    cell_ID_death = zeros(1, length(Lineage_2));%[cell_ID_death_daughter cell_ID_death(ind_empty)];
    cell_ID_birth = [cell_ID_birth_daughter cell_ID_birth(ind_empty)];
    cell_ID_plot = [cell_ID_plot_daughter cell_ID_plot(ind_empty)];
    if mod(i, 5) == 0
        disp(strcat("Iteration ", num2str(i), "/", num2str(Stat_t - 2)))
    end   
end
Table_fin = array2table(Tab_fin, 'VariableNames', {'Mother_ID', 'Daughter_ID_1', 'Daughter_ID_2', 'Lineage', 'Timepoint', 'Centroid_x', 'Centroid_y', 'orientation', 'axis_minor_length', 'axis_major_length', 'KeptOrNot', 'generation_time', 'generation', 'ID_death', 'ID_birth', 'GFP_cell_1', 'MCHE_cell_1', 'cell_ID_plot'});
Table_fin_tot = Table_fin;
Table_fin = Table_fin_tot;

%% Creation of a structure for each lineage

nb_Lineages = unique(Table_fin.Lineage);
nb_Lin = length(nb_Lineages);
Lineage_v2 = {};
Table_sorted = [Table_fin.Mother_ID Table_fin.Centroid_x Table_fin.Centroid_y Table_fin.Timepoint Table_fin.Lineage Table_fin.Daughter_ID_1 Table_fin.Daughter_ID_2 Table_fin.GFP_cell_1 Table_fin.MCHE_cell_1 Table_fin.cell_ID_plot Table_fin.generation];
Table_sorted = array2table(Table_sorted, 'VariableNames', {'MotherID',  'GeomX', 'GeomY', 'Timepoint', 'Lineage', 'DaughterID1', 'DaughterID2', 'GFP_mean', 'MCHE_mean', 'cell_ID_plot', 'generation'});%, 'orientation', 'axis_minor_length', 'axis_major_length', 'KeptOrNot', 'Life_time'});
Rows_NaN = Table_sorted.DaughterID2 == 0;
Table_sorted{Rows_NaN, 7} = NaN;

for i = 1:nb_Lin
    table_lineage_fin = Table_sorted(Table_sorted.Lineage == i, :);
    Lineage_v2{i} = table_lineage_fin;    
end

%% Creation gen_table

cellcounts_v1 = {};
[m, n] = size(Table_fin);
Table_fin.mu_max_vect = zeros(m, 1);

for i = 1:nb_Lin
    ind_temp_tot = find(Table_fin.Lineage == i);
    table_temp = Table_fin(Table_fin.Lineage == i, :);
    [m, n] = size(table_temp);
    table_temp.mu_max_vect = zeros(m, 1);
    [unique_ID_birth, ia] = unique(table_temp.ID_birth); [~, sorted_ia] = sort(ia); %The first apparition is sorted by timepoints
    unique_ID_birth = unique_ID_birth(sorted_ia);
    table_lineage_fin = [];
    [table_lineage_fin_tot, table_GFP_fin_tot, table_MCHE_fin_tot, table_plot_fin_tot, table_generation_fin_tot] = deal(nan(length(unique_ID_birth), Stat_t));
    for j = 1:length(unique_ID_birth)
        [temp_seq, temp_GFP, temp_MCHE, temp_plot, temp_generation] = deal([]);
        ind_temp = find(table_temp.ID_birth == unique_ID_birth(j));
        ind_find_mother = find(table_temp.Daughter_ID_1 == unique_ID_birth(j) | table_temp.Daughter_ID_2 == unique_ID_birth(j));
        temp = table_temp(ind_temp, :);
        temp = temp(table_temp.Timepoint(ind_temp) == max(table_temp.Timepoint(ind_temp)), :);
        mu_max_cell = log(2)/(1/3*length(ind_temp));
        table_temp.mu_max_vect(ind_temp) = mu_max_cell;
        Table_fin.mu_max_vect(ind_temp_tot(ind_temp)) = mu_max_cell;
        if ~isempty(ind_find_mother)
            mother_ind = table_temp.Mother_ID(ind_find_mother);
            last_time_point = (table_temp.Timepoint(ind_find_mother) + 1) - table_temp.Timepoint(1);
            previous_seq = find(table_lineage_fin_tot(:, last_time_point) == mother_ind); %Issue here
            previous_seq = previous_seq(1);
            temp_seq = table_lineage_fin_tot(previous_seq, :); temp_GFP = table_GFP_fin_tot(previous_seq, :); temp_MCHE = table_MCHE_fin_tot(previous_seq, :);
            temp_plot = table_plot_fin_tot(previous_seq, :); temp_generation = table_generation_fin_tot(previous_seq, :);
            temp_seq(isnan(temp_seq)) = []; temp_GFP(isnan(temp_plot)) = []; temp_MCHE(isnan(temp_plot)) = [];
            temp_plot(isnan(temp_plot)) = []; temp_generation(isnan(temp_generation)) = [];
        end
        table_lineage_fin = [table_lineage_fin; temp];
        temp_seq = [temp_seq table_temp.Mother_ID(ind_temp)']; temp_GFP = [temp_GFP table_temp.GFP_cell_1(ind_temp)']; temp_MCHE = [temp_MCHE table_temp.MCHE_cell_1(ind_temp)'];
        temp_plot = [temp_plot (table_temp.cell_ID_plot(ind_temp))']; temp_generation = [temp_generation (table_temp.generation(ind_temp) + 1)'];
        table_lineage_fin_tot(j, 1:length(temp_seq)) = temp_seq;
        table_GFP_fin_tot(j, 1:length(temp_seq)) = temp_GFP;
        table_MCHE_fin_tot(j, 1:length(temp_seq)) = temp_MCHE;
        table_plot_fin_tot(j, 1:length(temp_seq)) = temp_plot;
        table_generation_fin_tot(j, 1:length(temp_seq)) = temp_generation;
    end
    ind_not_dead = find(table_lineage_fin.ID_death == 0); table_lineage_fin.generation_time(ind_not_dead) = 0;
    cell_birth_tp = table_lineage_fin.Timepoint - table_lineage_fin.generation_time;
    cell_death_tp = table_lineage_fin.Timepoint;
    generation_time = table_lineage_fin.generation_time/3; %Every 20 minutes
    cell_death_tp(ind_not_dead) = nan; table_lineage_fin.ID_death(ind_not_dead) = nan; generation_time(ind_not_dead) = nan;
    table_lineage_fin = [table_lineage_fin.generation, table_lineage_fin.ID_birth, table_lineage_fin.ID_death, cell_birth_tp, cell_death_tp, generation_time];
    table_lineage_fin = array2table(table_lineage_fin, 'VariableNames', {'generation',  'ID_birth', 'ID_death', 'cell_birth_tp', 'cell_death_tp', 'generation_time'});
    cellcounts_v1.lineage(i).gen_table = table_lineage_fin;
    cellcounts_v1.lineage(i).mother_ID_table = table_lineage_fin_tot;
    cellcounts_v1.lineage(i).GFP_table = table_GFP_fin_tot;
    cellcounts_v1.lineage(i).MCHE_table = table_MCHE_fin_tot;
    cellcounts_v1.lineage(i).plot_table = table_plot_fin_tot;
    cellcounts_v1.lineage(i).generation_table = table_generation_fin_tot;
end


%% Part four: creation of the stationary phase tables
features = Data_kept;
data = sortrows(features,'Timepoint');
times = unique(data(:,4));
times = table2array(times);
Time_fin = max(unique(Table_fin.Timepoint)) - 0*4;
Fin_ind = find(Table_fin.Timepoint == Time_fin);

% go through the individual lineages
for i = 1:nb_Lin
    % define the stat phase mask of the cells    
    stat_mask = Lineage_v2{1,i}(Lineage_v2{1,i}.Timepoint == Time_fin, :); %Table_fin.Mother_ID(Fin_ind); %
    stat_tp = Time_fin; %max(Lineages{1,i}.Timepoint);
    
    %make new empty table with length of data to fill the IDs of next of kin
    %collect the data at each time point to make a new table that includes new cellIDs 
    stat_ID_table = stat_mask.MotherID; 
    stat_GFP_table = stat_mask.GFP_mean;
    stat_MCHE_table = stat_mask.MCHE_mean;
    stat_plot_table = stat_mask.cell_ID_plot;
    stat_generation_table = stat_mask.generation;
    
    stat_mask_geo = [stat_mask.GeomX, stat_mask.GeomY]; 
    
    %%go through all stationary phase time points and collect successively time=1, time+1 for distance comparison
    for j = stat_tp:(size(times, 1) - 2)
        table_lineage_fin = features(features.Timepoint == times(j+1), :);
        t = table2array(table_lineage_fin);%; data(data(:, 4) == times(j+1), :);
        
        %Modify the distance function
        Dist_temp = distEuclid([t(:,2),t(:,3)]', stat_mask_geo');
        [val, Ind_temp] = min(Dist_temp, [], 1);
            
        %columns in I are the IDs of the cells in stat_mask
        %the ones that are too far away become nan
    
        table_lineage_fin = table_lineage_fin(Ind_temp, :); %Table temp_reorganized
        ind_temp = find(val > 25);
        table_lineage_fin{ind_temp, :} = nan;
        Cell_IDs_temp = table2array(table_lineage_fin(:, 1));
        GFP_temp = table_lineage_fin.GFP_mean;
        MCHE_temp = table_lineage_fin.mCherry_mean;
        plot_temp = stat_mask.cell_ID_plot;
        generation_temp = stat_mask.generation;
        
        % identify the new cell IDs- can't use logic indexing because sometimes nan
        stat_ID_table = [stat_ID_table Cell_IDs_temp];
        stat_GFP_table = [stat_GFP_table GFP_temp];
        stat_MCHE_table = [stat_MCHE_table MCHE_temp];
        stat_plot_table = [stat_plot_table plot_temp];
        stat_generation_table = [stat_generation_table generation_temp];
    end
    
    %save into the cellcounts structure
    cellcounts_v1.lineage(i).stat_ID_table = stat_ID_table;
    cellcounts_v1.lineage(i).stat_GFP_table = stat_GFP_table;
    cellcounts_v1.lineage(i).stat_MCHE_table = stat_MCHE_table;
    cellcounts_v1.lineage(i).stat_plot_table = stat_plot_table;
    cellcounts_v1.lineage(i).stat_generation_table = stat_generation_table;
end

%% Connection to stationary phase

nb_Lin = length(nb_Lineages);
%Here is the issue with the stat time point. Check it
for i = 1:nb_Lin
    table_temp = Table_fin(Table_fin.Lineage == i, :); 
    mother_ID_temp = cellcounts_v1.lineage(i).mother_ID_table; GFP_temp = cellcounts_v1.lineage(i).GFP_table; MCHE_temp = cellcounts_v1.lineage(i).MCHE_table; 
    plot_temp = cellcounts_v1.lineage(i).plot_table; generation_temp = cellcounts_v1.lineage(i).generation_table;
    mother_ID_temp = mother_ID_temp(:, 1: stat_tp + 1); GFP_temp = GFP_temp(:, 1: stat_tp + 1); MCHE_temp = MCHE_temp(:, 1: stat_tp + 1);
    plot_temp = plot_temp(:, 1: stat_tp + 1); generation_temp = generation_temp(:, 1: stat_tp + 1);
    temp_stat = cellcounts_v1.lineage(i).stat_ID_table; GFP_temp_stat = cellcounts_v1.lineage(i).stat_GFP_table; MCHE_temp_stat = cellcounts_v1.lineage(i).stat_MCHE_table;
    plot_temp_stat = cellcounts_v1.lineage(i).stat_plot_table; generation_temp_stat = cellcounts_v1.lineage(i).stat_generation_table;
    if ~isempty(temp_stat)
        [n, m] = size(mother_ID_temp); [n_2, m_2] = size(temp_stat);
        [new_table, GFP_new_table, MCHE_new_table, plot_new_table, generation_new_table] = deal(nan(n, m + m_2 - 2));

        fin_ind = mother_ID_temp(:, stat_tp + 1 - table_temp.Timepoint(1));
        [val_temp, ind_1] = ismember(fin_ind, temp_stat(:, 1));
        ind_1 = ind_1(ind_1 > 0);

        cellcounts_v1 = Comb_stat_fun(cellcounts_v1, new_table, table_temp, mother_ID_temp, stat_tp, 'mother_ID_table', i, ind_1, val_temp, temp_stat);
        cellcounts_v1 = Comb_stat_fun(cellcounts_v1, GFP_new_table, table_temp, GFP_temp, stat_tp, 'GFP_table', i, ind_1, val_temp, GFP_temp_stat);
        cellcounts_v1 = Comb_stat_fun(cellcounts_v1, MCHE_new_table, table_temp, MCHE_temp, stat_tp, 'MCHE_table', i, ind_1, val_temp, MCHE_temp_stat);
        cellcounts_v1 = Comb_stat_fun(cellcounts_v1, plot_new_table, table_temp, plot_temp, stat_tp, 'plot_table', i, ind_1, val_temp, plot_temp_stat);
        cellcounts_v1 = Comb_stat_fun(cellcounts_v1, generation_new_table, table_temp, generation_temp, stat_tp, 'generation_table', i, ind_1, val_temp, generation_temp_stat);
    end
end

%% Remove lineages that are too small (optional)

for i = 1:nb_Lin
	if size(Lineage_v2{i}, 1) < 50
	    Lineage_v2{i} = [];
	    cellcounts_v1.lineage(i).gen_table = [];
	    cellcounts_v1.lineage(i).stat_ID_table = [];
	end
end

save(strcat('Data/', 'cellcount_v1.mat'), "cellcounts_v1")
save(strcat('Data/', 'Lineage_v2.mat'), "Lineage_v2")

%% example plot of the lineage tree
close all
k = 1;
Timepoint_display_fluo = 50; %Can be modified
lin_to_show = [1, 4, 5]; %1:nb_Lin %Insert the lineages to show. Warning: if there is a lot a lineage, don't show all of them, this will produce a lot of figures and Matlab can crash
for i = 1:length(lin_to_show)
    z = lin_to_show(i); %chose the lineage to plot    
    if ~isempty(Lineage_v2{i})
        figure(k)
        plot(1:size(cellcounts_v1.lineage(z).plot_table, 2), cellcounts_v1.lineage(z).plot_table, 'k');
        k = k + 1;

        figure(k)
        plot(1:size(cellcounts_v1.lineage(z).GFP_table,2), cellcounts_v1.lineage(z).GFP_table)
        k = k + 1;
        
        figure(k) %qqplot at time 50
        qqplot(cellcounts_v1.lineage(z).GFP_table(:, Timepoint_display_fluo))
        k = k + 1;

        figure(k)
        plot(1:size(cellcounts_v1.lineage(z).MCHE_table,2), cellcounts_v1.lineage(z).MCHE_table)
        k = k + 1;

        figure(k) %qqplot at time 50  
        qqplot(cellcounts_v1.lineage(z).MCHE_table(:, Timepoint_display_fluo))
        k = k + 1;
    end
end

%% Movie generation with cell shapes

close all;
real_sim = 1;
initial_time_point = 1; %0 if only stationary frame display, 1 if video generated
Data_file = Table_fin; %Name of the file and the name given to the video file if created
Data_file = Data_file(Data_file.KeptOrNot < 3, :);
name_movie_file = num_sheet;
nb_Lineages = unique(Table_fin.Lineage);
nb_Lin = length(nb_Lineages);
Tot_res = Fun_MovieSDEVec_Lineage(Data_file, real_sim, save_data, name_movie_file, initial_time_point);

initial_time_point = 1;% 0 if only stationary frame display, 1 if video generated
lin_to_show = 1:nb_Lin; %nb_Lin %Insert the lineages to show. 
for i = 1:length(lin_to_show)
    lineage_index = lin_to_show(i);
    % Movie generation or time point plot generation
    % Data_temp = Data_file(Data_file.Lineage == 5 | Data_file.Lineage == 4, :);
    Data_temp = Data_file(Data_file.Lineage == lineage_index, :);
    Tot_res = Fun_MovieSDEVec_Lineage(Data_temp, real_sim, save_data, name_movie_file, initial_time_point);
end

saveas(figure(1),strcat('Figures/', num_sheet,'_lineage_plot_.pdf'),'pdf')
writetable(Table_fin, strcat('Data/',num_sheet, '.xlsx'))