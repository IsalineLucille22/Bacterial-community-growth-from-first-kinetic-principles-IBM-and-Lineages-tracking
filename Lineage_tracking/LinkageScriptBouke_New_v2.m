%%%this script determines the best lineages that form microcolonies in a
%%%time-lapse experiment. It takes the output from Dimalis ('final_merged_table.xlsx') 
%%%and runs a set
%%%of calculations to find overlap between cells on consecutive images,
%%%from which then to determine the lineage. Lineages are only calculated
%%%based on the exponential phase. The last frame of the exponential phase
%%%is used as a mask to impose all stationary phase time points and connect
%%% those. The end result is a set of
%%%tables with cellIDs coupled to fluorescent data, organised per lineage
%%%and over time. In addition, the actual generation in the lineage is
%%%kept.
%%%Depending on the image quality and movement of individual images during
%%%the time-lapse, a correction can be attempted.
%%%An updated copy of the final_merged_table array will be saved that has
%%%all relevant fluorescence and cell measurements in the proper lineage
%%%order. This is specific for Bouke's data.

clear
close all


% %add the path where the scripts are
% addpath('/Users/jvanderm/Library/CloudStorage/OneDrive-UniversitédeLausanne/Lineage tracking/')
% 
% %Load the relevant data by specifying the path here
% cd('/Volumes/RECHERCHE/FAC/FBM/DMF/jvanderm/sinergia/D2c/RawData/Bouke/Microscopy/WORK_UNIL/PTYG_correct/20230321_Curtobacterium_PTYG_monoculture/Pos8')

% merged = readtable('final_merged_table.xlsx');
merged = readtable(strcat('Data/','final_merged_table_Bouke.xlsx'), 'Sheet', 1,'Format','auto'); 

Data_kept_init=[merged.Mask_nb, merged.Centroid_x, merged.Centroid_y, merged.Timepoint, merged.area, merged.eccentricity, merged.extent, merged.axis_minor_length, merged.axis_major_length, merged.feret_diameter_max, merged.orientation, merged.perimeter,merged.solidity, merged.Neighbors_FirstClosestDistance, merged.Neighbors_SecondClosestDistance, merged.GFP_mean, merged.GFP_median, merged.GFP_1st_quartile, merged.GFP_3rd_quartile, merged.GFP_sum, merged.GFP_min, merged.GFP_max,merged.mCherry_mean, merged.mCherry_median, merged.mCherry_1st_quartile, merged.mCherry_3rd_quartile, merged.mCherry_sum, merged.mCherry_min, merged.mCherry_max, merged.CFP_mean, merged.CFP_median, merged.CFP_1st_quartile, merged.CFP_3rd_quartile, merged.CFP_sum, merged.CFP_min, merged.CFP_max, merged.YFP_mean, merged.YFP_median, merged.YFP_1st_quartile, merged.YFP_3rd_quartile, merged.YFP_sum, merged.YFP_min, merged.YFP_max, merged.DAPI_mean, merged.DAPI_median, merged.DAPI_1st_quartile, merged.DAPI_3rd_quartile, merged.DAPI_sum, merged.DAPI_min, merged.DAPI_max];

Data_kept_init=array2table(Data_kept_init,'VariableNames',{'Mask_nb','Centroid_x','Centroid_y','Timepoint','area','eccentricity','extent','axis_minor_length','axis_major_length','feret_diameter_max','orientation','perimeter','solidity','Neighbors_FirstClosestDistance','Neighbors_SecondClosestDistance','GFP_mean','GFP_median','GFP_1st_quartile','GFP_3rd_quartile','GFP_sum','GFP_min','GFP_max','mCherry_mean','mCherry_median','mCherry_1st_quartile','mCherry_3rd_quartile','mCherry_sum','mCherry_min','mCherry_max','CFP_mean','CFP_median','CFP_1st_quartile','CFP_3rd_quartile','CFP_sum','CFP_min','CFP_max','YFP_mean','YFP_median','YFP_1st_quartile','YFP_3rd_quartile','YFP_sum','YFP_min','YFP_max','DAPI_mean','DAPI_median','DAPI_1st_quartile','DAPI_3rd_quartile','DAPI_sum','DAPI_min','DAPI_max'});
num_sheet = 'Pos8';

Data_kept_init.Mask_nb = (1:height(Data_kept_init))'; %reinitialize cell numbers
Data_kept_init.Timepoint = Data_kept_init.Timepoint - min(Data_kept_init.Timepoint); %Modify time to start at 0 anyway
real_sim = 1; %1 if real, 0 if sim
[~, sorted_time] = sort(Data_kept_init.Timepoint);
Data_kept_init = Data_kept_init(sorted_time, :);
Data_kept_init.orientation = -(pi/2 - Data_kept_init.orientation);
Data_kept_init.Centroid_y = 2050 - Data_kept_init.Centroid_y;

%Sometimes, there are missing time points, which can cause issues later. To address this, we force time points to increment only by 1.
Data_kept_init = Arrage_Timepoint_fun(Data_kept_init); 
time_point = unique(Data_kept_init.Timepoint);
time_point_init = 1;
Stat_t = Stat_Phase_fun(Data_kept_init);

%if a correction for moving images needs to be made, keep the following
%line. Otherwise inactivate it.
Data_kept_init = fun_Slinding_Effect_v2(Data_kept_init);
weights = [0.25, 0.25, 0.25, 0.25]; %Not used

%if your data require a more elaborate overlap calculation set the
%following line to 0; otherwise chose 1. More elaborate calculation takes
%a bit more time!

method_number = 0;
if method_number == 0
    real_sim = -1; %1 if real, 0 if sim
end

len_time_point = Stat_t;

%Set the range to test increasing cell diameters to find potential overlap.
%Typically, a range of 1:12 is okay.

overlap_val_margin_tot = 1:12; 
nb_No_unique = inf*ones(1, length(overlap_val_margin_tot)); 
Table_fin_tot = {};
int_num_cells = length(find(Data_kept_init.Timepoint == min(time_point)));
new_column = ones(height(Data_kept_init), 1);  %Generate new data for the column
Data_kept_init.KeptOrNot = new_column;  %Assign the new column
new_column = zeros(height(Data_kept_init), 1);  %Generate new data for the column
Data_kept_init.Life_time = new_column;  %Assign the new column

mean_length =  mean(Data_kept_init.axis_major_length(~isnan(Data_kept_init.axis_major_length))); %mean(Data_1.axis_major_length(~isnan(Data_1.axis_major_length)));%
mean_d =  mean(Data_kept_init.axis_minor_length(~isnan(Data_kept_init.axis_minor_length)));%mean(Data_1.axis_minor_length(~isnan(Data_1.axis_minor_length)));%

% The following lines are intended to be able to correct for consistency of
% cell shape. Sometimes cells have aberrant shapes caused by the
% segmentation. We envisioned four types of correction. (i) remove objects that are too small
% and probably artefacts. (ii) Split cell objects into two when they are
% clearly too long. (iii) Correct the thickness of cells that are below or
% above the mean measured diameter to that mean, or (iv) for the length of
% cells.
%
% Data_kept_init = Remove_Cells(Data_kept_init, 'axis_minor_length', mean_d, 0.1, 1.2); %case (i)
% Data_kept_init = Length_cells_mod(Data_kept_init, mean_length, 4); %case(ii) 
% Data_kept_init = Correction_Cell_Sizes(Data_kept_init, 'axis_minor_length', mean_d, 0.8, 1.4); %case (iii)
% Data_kept_init = Correction_Cell_Sizes(Data_kept_init, 'axis_major_length', mean_length, 0.4, 2.5); % case (iv)
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

%%the program will now loop over all time points and will determine both at
%%the next and the after-next time point whether overlap exists. It will
%%also attempt to replace a cell should such object be missing from a next
%%time point.
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
    [Data_kept, Pos_X_1] = Fill_NaN_Pos(ind_t_1, Data_kept, Pos_X_1);
    Pos_X_2 = [Data_2.Centroid_x, Data_2.Centroid_y]';
    [Data_kept, Pos_X_2] = Fill_NaN_Pos(ind_t_2, Data_kept, Pos_X_2);
    Pos_X_3 = [Data_3.Centroid_x, Data_3.Centroid_y]';
    [Data_kept, Pos_X_3] = Fill_NaN_Pos(ind_t_3, Data_kept, Pos_X_3);
    Dist = distEuclid(Pos_X_1, Pos_X_2);
    Dist_23 = distEuclid(Pos_X_2, Pos_X_3);
    Dist_23 = [Dist' Dist_23];

    vect_Cell_length_1 = Data_1.axis_major_length; 
    height_cell_1 = Data_1.axis_minor_length;
    cell_angle_1 = Data_1.orientation;
    Kept_NOR_1 = Data_1.KeptOrNot;
    cell_1_area = Data_1.area; cell_1_eccentricity = Data_1.eccentricity; cell_1_extent = Data_1.extent; cell_1_feret_diameter_max = Data_1.feret_diameter_max; cell_1_perimeter = Data_1.perimeter; cell_1_solidity = Data_1.solidity; cell_1_GFP_mean = Data_1. GFP_mean; cell_1_GFP_median = Data_1.GFP_median; cell_1_GFP_1st_quartile = Data_1.GFP_1st_quartile; cell_1_GFP_3rd_quartile = Data_1.GFP_3rd_quartile; cell_1_GFP_sum = Data_1.GFP_sum; cell_1_GFP_min = Data_1.GFP_min; cell_1_GFP_max = Data_1.GFP_max; cell_1_mCherry_mean = Data_1. mCherry_mean; cell_1_mCherry_median = Data_1.mCherry_median; cell_1_mCherry_1st_quartile = Data_1. mCherry_1st_quartile; cell_1_mCherry_3rd_quartile = Data_1.mCherry_3rd_quartile; cell_1_mCherry_sum = Data_1.mCherry_sum; cell_1_mCherry_min = Data_1.mCherry_min; cell_1_mCherry_max = Data_1.mCherry_max; cell_1_CFP_mean = Data_1.CFP_mean; cell_1_CFP_median = Data_1.CFP_median; cell_1_CFP_1st_quartile = Data_1.CFP_1st_quartile; cell_1_CFP_3rd_quartile = Data_1.CFP_3rd_quartile; cell_1_CFP_sum = Data_1.CFP_sum; cell_1_CFP_min = Data_1.CFP_min; cell_1_CFP_max = Data_1.CFP_max; cell_1_YFP_mean = Data_1.YFP_mean; cell_1_YFP_median = Data_1. YFP_median; cell_1_YFP_1st_quartile = Data_1.YFP_1st_quartile; cell_1_YFP_3rd_quartile = Data_1.YFP_3rd_quartile; cell_1_YFP_sum = Data_1.YFP_sum; cell_1_YFP_min = Data_1.YFP_min; cell_1_YFP_max = Data_1.YFP_max; cell_1_DAPI_mean = Data_1.DAPI_mean; cell_1_DAPI_median = Data_1.DAPI_median; cell_1_DAPI_1st_quartile = Data_1.DAPI_1st_quartile; cell_1_DAPI_3rd_quartile = Data_1.DAPI_3rd_quartile; cell_1_DAPI_sum = Data_1.DAPI_sum; cell_1_DAPI_min = Data_1.DAPI_min; cell_1_DAPI_max=Data_1.DAPI_max;
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

    % here is reference to the overlap function method

    h_mat = H_mat_fun(Dist, n_1, n_2, vect_Cell_length_1 - real_sim*height_cell_1, vect_Cell_length_2 - real_sim*height_cell_2, height_cell_1, height_cell_2, Seg_tot_1, Seg_tot_2, weights, Seg_tot_1, Seg_tot_2, method_number); 
    h_mat(h_mat < 0) = 0;
    [val, Lineage_ind] = max(h_mat, [], 1);
    %%%
    h_mat_1 = ones(n_1, n_2); h_test = h_mat;
    while m < length(overlap_val_margin_tot) && sum(val_temp == 0) > 0
    
        %%%
        col_inf = find(isinf(sum(h_mat_1)));
        col_rem = find(sum(h_mat_1) < n_1);
        while isequal((val == 0), (val_temp == 0)) && m < length(overlap_val_margin_tot)
            overlap_val_margin = overlap_val_margin_tot(m);
            h_mat = H_mat_fun(Dist, n_1, n_2, vect_Cell_length_1 - real_sim*height_cell_1, vect_Cell_length_2 - real_sim*height_cell_2, overlap_val_margin*height_cell_1, overlap_val_margin*height_cell_2, Seg_tot_1, Seg_tot_2, weights, Seg_tot_1_rect, Seg_tot_2_rect, method_number); 
            h_mat(h_mat < 0) = 0;
            h_test = h_mat;
            h_mat(:, col_inf) = h_mat_1(:, col_inf);
            h_mat(:, col_rem) = h_mat_1(:, col_rem).*h_test(:, col_rem);
            h_test(:, col_rem) = h_mat_1(:, col_rem).*h_test(:, col_rem);
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
                        ind_temp = ind_temp(b);
                        h_mat_1(j, Ind_Daughters_NK) = 0; h_test(j, Ind_Daughters_NK) = 0;
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
                if ~isempty(ID_Daughters)
                    Daughter_temp(j, :) = ID_Daughters;
                end
                %%%
                h_mat_1(1:n_1~=j, ind_temp) = 0; h_mat_1(j, ind_temp) = inf;
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
        m = m + 1;
        val_temp = val;
    end
    inf_h_mat = find(isinf(h_mat));
    h_test(inf_h_mat) = inf;
    h_mat = h_test;

    %%%
    Daughter_temp_1 = Daughter_temp;
    Daughter_temp_1(Daughter_temp(:, 1) - Daughter_temp(:, 2) == 0, 2) = 0;
    Daughters_vect = [Daughter_temp_1(:,1); Daughter_temp_1(:,2)];
    Daughters_vect = Daughters_vect(Daughters_vect ~= 0);
    unique_vals_Daughters = unique(Daughters_vect);
    %%%

    %%%% Calculate the life time for each cell
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
    Lineage_2(val == 0) = max_Lineage + 1:max_Lineage + nb_zeros;
    
    % Add missing/disappearing cells during one time frame
    [Data_kept, Daughter_temp, ind_empty] = Disappearing_Cells(Data_kept, Daughter_temp, ind_empty, ind_t_1, t_2);
    
    % Calculate the actual generation time from the time of birth and death
    % of each cell in the lineage

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
    
    Daughter_temp(Daughter_temp(:, 1) - Daughter_temp(:, 2) == 0, 2) = 0;
    Tab_1 = [cell_ID_1, Daughter_temp, Lineage_1', Data_1.Timepoint,  Pos_X_1', cell_angle_1, height_cell_1, vect_Cell_length_1, Kept_NOR_1, Life_time', generation_time_save', cell_ID_death', cell_ID_birth', cell_ID_plot', cell_1_area, cell_1_eccentricity, cell_1_extent, cell_1_feret_diameter_max, cell_1_perimeter, cell_1_solidity, cell_1_GFP_mean, cell_1_GFP_median, cell_1_GFP_1st_quartile, cell_1_GFP_3rd_quartile, cell_1_GFP_sum, cell_1_GFP_min, cell_1_GFP_max, cell_1_mCherry_mean, cell_1_mCherry_median, cell_1_mCherry_1st_quartile, cell_1_mCherry_3rd_quartile, cell_1_mCherry_sum, cell_1_mCherry_min, cell_1_mCherry_max, cell_1_CFP_mean, cell_1_CFP_median, cell_1_CFP_1st_quartile, cell_1_CFP_3rd_quartile, cell_1_CFP_sum, cell_1_CFP_min, cell_1_CFP_max, cell_1_YFP_mean, cell_1_YFP_median, cell_1_YFP_1st_quartile, cell_1_YFP_3rd_quartile, cell_1_YFP_sum, cell_1_YFP_min, cell_1_YFP_max, cell_1_DAPI_mean, cell_1_DAPI_median, cell_1_DAPI_1st_quartile, cell_1_DAPI_3rd_quartile, cell_1_DAPI_sum, cell_1_DAPI_min, cell_1_DAPI_max]; %Create here the all_tracks table 
    Tab_fin = [Tab_fin; Tab_1];
    max_Lineage = max_Lineage + nb_zeros;
    [~, sorted_time] = sort(Data_kept.Timepoint);
    Data_kept = Data_kept(sorted_time, :);
    Lineage_2 = [Lineage_2 Lineage_1(ind_empty)];
    Life_time_daughter = [Life_time_daughter Life_time(ind_empty)];
    
    generation_time = [generation_time_daughter generation_time(ind_empty)];
    cell_ID_death = zeros(1, length(Lineage_2));%[cell_ID_death_daughter cell_ID_death(ind_empty)];
    cell_ID_birth = [cell_ID_birth_daughter cell_ID_birth(ind_empty)];
    cell_ID_plot = [cell_ID_plot_daughter cell_ID_plot(ind_empty)];
    if mod(i, 5) == 0
        disp(strcat("Iteration ", num2str(i), "/", num2str(Stat_t - 2)))
    end   
end
Table_fin = array2table(Tab_fin, 'VariableNames', {'Mother_ID', 'Daughter_ID_1', 'Daughter_ID_2', 'Lineage', 'Timepoint', 'Centroid_x', 'Centroid_y', 'orientation', 'axis_minor_length', 'axis_major_length', 'KeptOrNot', 'generation_time', 'generation', 'ID_death', 'ID_birth','cell_ID_plot','area','eccentricity','extent','feret_diameter_max','perimeter','solidity','GFP_mean','GFP_median','GFP_1st_quartile','GFP_3rd_quartile','GFP_sum','GFP_min','GFP_max','mCherry_mean','mCherry_median','mCherry_1st_quartile','mCherry_3rd_quartile','mCherry_sum','mCherry_min','mCherry_max','CFP_mean','CFP_median','CFP_1st_quartile','CFP_3rd_quartile','CFP_sum','CFP_min','CFP_max','YFP_mean','YFP_median','YFP_1st_quartile','YFP_3rd_quartile','YFP_sum','YFP_min','YFP_max','DAPI_mean','DAPI_median','DAPI_1st_quartile','DAPI_3rd_quartile','DAPI_sum','DAPI_min','DAPI_max'});
Table_fin_tot = Table_fin;
Table_fin = Table_fin_tot;

%% Creation of a structure for each lineage

nb_Lineages = unique(Table_fin.Lineage);
nb_Lin = length(nb_Lineages);
Lineage_v2 = {};
Table_sorted = [Table_fin.Mother_ID Table_fin.Centroid_x Table_fin.Centroid_y Table_fin.Timepoint Table_fin.Lineage Table_fin.Daughter_ID_1 Table_fin.Daughter_ID_2 Table_fin.GFP_mean Table_fin.mCherry_mean Table_fin.cell_ID_plot Table_fin.generation];
Table_sorted = array2table(Table_sorted, 'VariableNames', {'MotherID',  'GeomX', 'GeomY', 'Timepoint', 'Lineage', 'DaughterID1', 'DaughterID2', 'GFP_mean', 'MCHE_mean', 'cell_ID_plot', 'generation'});%, 'orientation', 'axis_minor_length', 'axis_major_length', 'KeptOrNot', 'Life_time'});
Rows_NaN = Table_sorted.DaughterID2 == 0;
Table_sorted{Rows_NaN, 7} = NaN;

for i = 1:nb_Lin
    table_lineage_fin = Table_sorted(Table_sorted.Lineage == i, :);
    Lineage_v2{i} = table_lineage_fin;    
end

%% Creation of the cellcounts structure
% that collects the generation table, the cell ID table, the plot table
% and the fluorescence tables

cellcounts_v1 = {};

%Addition of the growth rate in the Table_fin
[m, n] = size(Table_fin);
Table_fin.mu_max_vect = zeros(m, 1);

for i = 1:nb_Lin
    ind_temp_tot = find(Table_fin.Lineage == i);
    table_temp = Table_fin(Table_fin.Lineage == i, :);  
    [unique_ID_birth, ia] = unique(table_temp.ID_birth); [~, sorted_ia] = sort(ia); %The first apparition is sorted by timepoints
    unique_ID_birth = unique_ID_birth(sorted_ia);
    table_lineage_fin = [];
    [table_lineage_fin_tot, table_GFP_fin_tot, table_MCHE_fin_tot, table_plot_fin_tot, table_generation_fin_tot] = deal(nan(length(unique_ID_birth), Stat_t));
    for j = 1:length(unique_ID_birth)
        [temp_seq, temp_GFP, temp_MCHE, temp_plot, temp_generation] = deal([]);
        % temp_seq = []; temp_GFP = []; temp_plot = []; temp_generation = [];
        ind_temp = find(table_temp.ID_birth == unique_ID_birth(j));
        ind_find_mother = find(table_temp.Daughter_ID_1 == unique_ID_birth(j) | table_temp.Daughter_ID_2 == unique_ID_birth(j));
        temp = table_temp(ind_temp, :);
        temp = temp(table_temp.Timepoint(ind_temp) == max(table_temp.Timepoint(ind_temp)), :);
        mu_max_cell = log(2)/(1/6*length(ind_temp));
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
        temp_seq = [temp_seq table_temp.Mother_ID(ind_temp)']; temp_GFP = [temp_GFP table_temp.GFP_mean(ind_temp)']; temp_MCHE = [temp_MCHE table_temp.mCherry_mean(ind_temp)'];
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
    generation_time = table_lineage_fin.generation_time/6; %Every 10 minutes
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


%% Creation of the stationary phase tables

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
    if ~isempty(stat_mask_geo)
    
    %%go through all stationary phase time points and collect successively time=1, time+1 for distance comparison
        for j = stat_tp:(size(times, 1) - 2)
            table_lineage_fin = features(features.Timepoint == times(j+1), :);
            t = table2array(table_lineage_fin);%; data(data(:, 4) == times(j+1), :);
            
            %Modify the distance function
            Dist_temp = distEuclid([t(:,2),t(:,3)]', stat_mask_geo');
            [val, Ind_temp] = min(Dist_temp, [], 1);
    
            % Ind_temp(2) = Ind_temp(4);
            unique_val = unique(Ind_temp);
            counts = histcounts(Ind_temp, [unique_val, max(unique_val) + 1]);
            while max(counts) > 1
                ind_mult_val = find(counts > 1);
                val_temp = find(Ind_temp == unique_val(ind_mult_val(1)));
                [val_min, min_pos] = min(val(val_temp));
                max_pos = val_temp(setdiff(1:length(val_temp), min_pos));
                Dist_temp(unique_val(ind_mult_val(1)), max_pos) = inf;
                [val, Ind_temp] = min(Dist_temp, [], 1);
                unique_val = unique(Ind_temp);
                counts = histcounts(Ind_temp, [unique_val, max(unique_val) + 1]);
            end
                
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
end

%% Connect exponential phase to stationary phase data

nb_Lin = length(nb_Lineages);
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
	if size(Lineage_v2{i}, 1) < 10
	    Lineage_v2{i} = [];
	    cellcounts_v1.lineage(i).gen_table = [];
	    cellcounts_v1.lineage(i).stat_ID_table = [];
	end
end

save(strcat('Data/', num_sheet, 'cellcount_v1.mat'), "cellcounts_v1")
save(strcat('Data/', num_sheet, 'Lineage_v2.mat'), "Lineage_v2")

%% example plot of the lineage tree
close all
k = 1;
Timepoint_display_fluo = 40; %Can be modified
lin_to_show = 1:5; %1:nb_Lin %Insert the lineages to show. Warning: if there is a lot a lineage, don't show all of them, this will produce a lot of figures and Matlab can crash
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
%Save or Not
save_data = 1; %1 if the video is saved, 0 otherwise

close all;
real_sim = 1;
initial_time_point = 0;% 0 if only stationary frame display, 1 if video generated
Data_file = Table_fin;%Name of the file and the name given to the video file if created
Data_file = Data_file(Data_file.KeptOrNot < 2, :); %Remove cells that are manually kept x time. If inf kept all cells even the ones that are manually kept, if 1 display only cells detected by Dimalis
name_movie_file = num_sheet;
nb_Lineages = unique(Table_fin.Lineage);
nb_Lin = length(nb_Lineages);
% Tot_res = Fun_MovieSDEVec_Lineage(Data_file, real_sim, save_data, name_movie_file, initial_time_point); %Colored by lineages
% Tot_res_mu = Fun_MovieSDEVec_Growth_rate(Data_file, real_sim, save_data, name_movie_file, 1); %Colored by growth rates
Tot_res = Fun_SpecificTime_Growth_rates(Data_file, real_sim, save_data, 15);
% to include the lineage number on the image, inactivate the line starting
% with text(P_temp_S(1, :) in the function script.
% save the figure
saveas(figure(1),strcat('Figures/', num_sheet,'_lineage_plot_.pdf'),'pdf')
writetable(Table_fin, strcat('Data/','lineage_table_',num_sheet, '.xlsx'))


%% Individual lineage
initial_time_point = 0;% 0 if only stationary frame display, 1 if video generated
lin_to_show = 3; %1:nb_Lin; %nb_Lin %Insert the lineages to show. 
for i = 1:length(lin_to_show)
    lineage_index = lin_to_show(i);
    % Movie generation or time point plot generation
    Data_temp = Data_file(Data_file.Lineage == lineage_index, :);
    Tot_res = Fun_MovieSDEVec_Lineage(Data_temp, real_sim, save_data, name_movie_file, initial_time_point);
end