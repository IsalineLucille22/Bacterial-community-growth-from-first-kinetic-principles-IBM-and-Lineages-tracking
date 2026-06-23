clear
close all


% Define initial weights (e.g., equal weights as a starting point)
initial_weights = 0.*[1/3, 1/3, 1/3]; % Adjust based on the number of weights
weights(1)  = pi/4;

% % Define bounds for weights (e.g., [0, 1] if weights should be between 0 and 1)
% lb = [0, 0, 0]; % Lower bounds
% ub = [1, 1, 1]; % Upper bounds
% 
% % Use fmincon to optimize weights
% options = optimoptions('fmincon', 'MaxIterations', 1000, 'Display', 'iter');
% optimal_weights = fmincon(@(w) WeightsComp(w), initial_weights, [], [], [], [], lb, ub, [], options);
% 
% % Display results
% disp('Optimal Weights:');
% disp(optimal_weights);
cost = zeros(1, 20);
for i = 1:20
    cost(i) = WeightsComp(initial_weights);
    initial_weights = initial_weights + pi/10;
%     initial_weights = max(initial_weights + normrnd(0, 0.1, 1, 3), 0);
%     initial_weights = initial_weights/sum(initial_weights);
end


function cost = WeightsComp(weights)

%Loadind data
num_sheet = 6; 'Anthony_Pos1';
real_sim = 0; %1 if real, 0 if sim
Data_kept_init = readtable(strcat('Data/','Simulated_Data_v2.xlsx'), 'Sheet', num_sheet,'Format','auto'); 
Data_kept_init.Timepoint = Data_kept_init.Timepoint - min(Data_kept_init.Timepoint); %Modify time to start at 0 anyway
Data_kept_init = Arrage_Timepoint_fun(Data_kept_init);
% [Data_kept_init, translated_vector] = AdditionNoise(Data_kept_init);
[~, sorted_time] = sort(Data_kept_init.Timepoint);
Data_kept_init = Data_kept_init(sorted_time, :);
time_point = unique(Data_kept_init.Timepoint);
time_point_init = 1;
Stat_t = length(time_point);%40;%
% Data_kept_init = fun_Slinding_Effect_v2(Data_kept_init);

overlap_val_margin_tot = 1:6; %Sequence for the overlapping margings (intentional increase of the cell diameter)
nb_No_unique = inf*ones(1, length(overlap_val_margin_tot)); Table_fin_tot = {};
int_num_cells = length(find(Data_kept_init.Timepoint == min(time_point)));
new_column = ones(height(Data_kept_init), 1);  %Generate new data for the column
Data_kept_init.KeptOrNot = new_column;  %Assign the new column
new_column = zeros(height(Data_kept_init), 1);  %Generate new data for the column
Data_kept_init.Life_time = new_column;  %Assign the new column

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
for i = time_point_init:(Stat_t - 2)
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
    vect_Cell_length_2 = Data_2.axis_major_length; 
    height_cell_2 = Data_2.axis_minor_length;
    cell_angle_2 = Data_2.orientation;
    Kept_NOR_2 = Data_2.KeptOrNot;
    vect_Cell_length_3 = Data_3.axis_major_length; 
    height_cell_3 = Data_3.axis_minor_length;
    cell_angle_3 = Data_3.orientation;

    %Fill the NaN values
    %Combination t-1 and t+1 for t_2
    [Data_kept, vect_Cell_length_1, height_cell_1, cell_angle_1] = Fill_NaN(ind_t_1, Dist', Data_kept, vect_Cell_length_1, height_cell_1, cell_angle_1, Pos_X_1, vect_Cell_length_2, height_cell_2, cell_angle_2, []);
    [Data_kept, vect_Cell_length_2, height_cell_2, cell_angle_2] = Fill_NaN(ind_t_2, Dist_23', Data_kept, vect_Cell_length_2, height_cell_2, cell_angle_2, Pos_X_2, [vect_Cell_length_1; vect_Cell_length_3], [height_cell_1; height_cell_3], [cell_angle_1; cell_angle_3], []);

    vect_Cell_length_1 = max(vect_Cell_length_1, mean(Data_kept.axis_major_length(ind_t_1)));
    n_1 = nb_Cell_t_1; %height(Data_1(:, 1));
    Seg_tot_1 = arrayfun(@(x) Rect2Seg([Pos_X_1(1,x) Pos_X_1(2,x) (vect_Cell_length_1(x) - real_sim*height_cell_1(x)) height_cell_1(x)], cell_angle_1(x)),1:n_1,'UniformOutput',false); %Find an alternative
    cell_ID_1 = Data_1.Mask_nb;
    height_cell_1 = Data_kept.axis_minor_length(ind_t_1);%mean(Data_kept.axis_minor_length(ind_t_1))*ones(nb_Cell_t_1, 1);
    
    vect_Cell_length_2 = max(vect_Cell_length_2, mean(Data_kept.axis_major_length(ind_t_2)));
    n_2 = nb_Cell_t_2; %height(Data_2(:, 1)); 
    Seg_tot_2 = arrayfun(@(x) Rect2Seg([Pos_X_2(1,x) Pos_X_2(2,x) (vect_Cell_length_2(x) - real_sim*height_cell_2(x)) height_cell_2(x)], cell_angle_2(x)),1:n_2,'UniformOutput',false); %Find an alternative
    cell_ID_2 = Data_2.Mask_nb;
    height_cell_2 = Data_kept.axis_minor_length(ind_t_2);%mean(Data_kept.axis_minor_length(ind_t_2))*ones(nb_Cell_t_2, 1);
    
    m = 1; val_temp = 0;
    Lineage_1 = Lineage_2;
    Life_time = Life_time_daughter;
    Daughter_temp = zeros(n_1, 2);
    linear_indices = [];
    h_mat = H_mat_fun(Dist, n_1, n_2, vect_Cell_length_1 - real_sim*height_cell_1, vect_Cell_length_2 - real_sim*height_cell_2, height_cell_1, height_cell_2, Seg_tot_1, Seg_tot_2, weights); 
    h_mat(h_mat < 0) = 0;
    h_mat(linear_indices) = inf;
    [val, Lineage_ind] = max(h_mat, [], 1);
    while m < length(overlap_val_margin_tot) && sum(val_temp == 0) > 0
        overlap_val_margin = overlap_val_margin_tot(m);
    
        while isequal((val == 0), (val_temp == 0)) && m < length(overlap_val_margin_tot)
            h_mat = H_mat_fun(Dist, n_1, n_2, vect_Cell_length_1 - real_sim*height_cell_1, vect_Cell_length_2 - real_sim*height_cell_2, overlap_val_margin*height_cell_1, overlap_val_margin*height_cell_2, Seg_tot_1, Seg_tot_2, weights); 
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

        %%%%
        ind_empty = [];
        for j = 1:n_1
            ind_temp = find(Lineage_ind == j);
            ID_Daughters = cell_ID_2(ind_temp);
            if length(ID_Daughters) > 2
                [~, b] = maxk(h_mat(j, ind_temp), 2);
                ID_Daughters = ID_Daughters(b);
            end
            if ~isempty(ID_Daughters)
                Daughter_temp(j, :) = ID_Daughters;
            else 
                ind_empty = [ind_empty, j]; %For empty mother look at t + 1 and add them if they have daugthers 
            end
            [val, Lineage_ind] = max(h_mat, [], 1);
            Lineage_2 = Lineage_1(Lineage_ind);
            nb_zeros = sum(val == 0);
            Lineage_2(val == 0) = max_Lineage + 1:max_Lineage + nb_zeros;
            Lineage_ind(val == 0) = 0;
        end

        ind_daughter = 1:n_2;
        ind_daughter(val == 0) = [];
        if ~isempty(ind_daughter)
            linear_indices = sub2ind(size(h_mat), Lineage_ind(Lineage_ind ~= 0), ind_daughter);  %Replace by inf the indices corresponding to a relation mother-daughter attributed
        end
        m = m + 1;
        val_temp = val;
    end

    %%%% Life_time
    [h_mat, Daughter_temp, Life_time, Life_add, ind_empty] = Life_Time_fun_v2(h_mat, Daughter_temp, Life_time, Lineage_ind, cell_ID_2);
    Life_time_temp = Life_time;
    Life_time_temp(Life_add == 0) = 0;
    [val, Lineage_ind] = max(h_mat, [], 1);
    Life_time_daughter = Life_time_temp(Lineage_ind);
    Life_time_daughter(val == 0) = 0;
    Lineage_2 = Lineage_1(Lineage_ind);
    nb_zeros = sum(val == 0);
    Lineage_2(val == 0) = max_Lineage + 1:max_Lineage + nb_zeros;
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
    Tab_1 = [cell_ID_1, Daughter_temp, Lineage_1', Data_1.Timepoint,  Pos_X_1', cell_angle_1, height_cell_1, vect_Cell_length_1, Kept_NOR_1, Life_time', generation_time_save', cell_ID_death', cell_ID_birth', cell_ID_plot']; %Create here the all_tracks table 
    Tab_fin = [Tab_fin; Tab_1];
    max_Lineage = max_Lineage + nb_zeros;
    [~, sorted_time] = sort(Data_kept.Timepoint);
    Data_kept = Data_kept(sorted_time, :);
    generation_time = generation_time_daughter; %[generation_time_daughter generation_time(ind_empty)];%generation_time_daughter; %
    cell_ID_death = cell_ID_death_daughter; %[cell_ID_death_daughter cell_ID_death(ind_empty)];%zeros(1, length(Lineage_2));%
    cell_ID_birth = cell_ID_birth_daughter;%[cell_ID_birth_daughter cell_ID_birth(ind_empty)];%cell_ID_birth_daughter;%
    cell_ID_plot = cell_ID_plot_daughter;%[cell_ID_plot_daughter cell_ID_plot(ind_empty)];%cell_ID_plot_daughter;%
end
Table_fin = array2table(Tab_fin, 'VariableNames', {'Mother_ID', 'Daughter_ID_1', 'Daughter_ID_2', 'Lineage', 'Timepoint', 'Centroid_x', 'Centroid_y', 'orientation', 'axis_minor_length', 'axis_major_length', 'KeptOrNot', 'generation_time', 'generation', 'ID_death', 'ID_birth', 'cell_ID_plot'});
Table_fin_tot = Table_fin;
Table_fin = Table_fin_tot;

%%% Comparison lineage

Data_kept(Data_kept.Timepoint > t_1, :) = [];

cost = weight_cost(Data_kept, Table_fin, weights);
disp('Current weights:');
disp(weights);
disp(['Current cost: ', num2str(cost)]);
end