%this script takes the final biomass or area data from simulations to compare
% with the actual observed species biomass proportions on SE in the full minisyncom
% This variant uses simulations with mu's from the monocultures (files named 'SynCom7_SE')
% of with the actual mu's for each of the target cultures measured in the mini-SynCom (files named 'SynCom7_SE_mix_mu')
% Change the corresponding lines!

clear;
close all;


%%%%%Initialization cells positions based on real data%%%%%

path_init = '/Users/iguex/Library/CloudStorage/OneDrive-UniversitédeLausanne/Images_analysis/Scripts/Zenodo Paper/Spatial scripts/';
path_Data_position = '/Volumes/RECHERCHE/FAC/FBM/DMF/jvanderm/sinergia/D2c/RawData/Tania/Microscopy/SE_experiments/20251022_minisyncom_SE/Dimalis';
% mu_Mody = readtable(strcat('/Volumes/RECHERCHE/FAC/FBM/DMF/jvanderm/sinergia/D2c/UserData/Isaline/Model/', 'mean-observed-paired-slopes.csv')); %Slope in the 2nd column indicates the effect of the 1st sepcies on the growth rate of the 2nd species (so it is the max umax of the 2nd sepcies that is modified)
mu_Mody = readtable(strcat('/Volumes/RECHERCHE/FAC/FBM/DMF/jvanderm/sinergia/D2c/UserData/Isaline/Model/', 'yield_correction_table.csv')); %Yields modification
Data_path = strcat(path_init, 'Data/');
cd(path_Data_position)
init_name = 'Zenodo_SynCom7_SE_mono_growth_rate_with_umax_higher_general_yield_0.6';
addpath(path_init)
addpath(path_Data_position)
addpath(strcat(path_init, '/Data'))
color_set = [
    1 0 0;   %red
    0 0 1;   %blue
    1 0.5 0; %orange
    0 0 0;   %black
    0 1 0;   %green
    0 1 1;   %cyan
    1 0 1;   %magenta
];

positions = dir('Pos*');

fig2 = figure(2); clf(fig2)
ax2 = axes(fig2); hold(ax2,'on')
xlabel(ax2,'Simulated area')
ylabel(ax2,'Observed area')
xlim(ax2,[0 3e05])
ylim(ax2,[0 3e05])
fplot(ax2, @(x) x, [0 3e05], 'k--', 'Linewidth', 0.5)
fig3 = figure(3); clf(fig3)
ax3 = axes(fig3); hold(ax3,'on')
xlabel(ax3,'Simulated area')
ylabel(ax3,'Observed area')
xlim(ax3,[0 1e03])
ylim(ax3,[0 1e03])
fplot(ax3, @(x) x, [0 1e03], 'k--', 'Linewidth', 0.5)
nb_pos = length(positions);
[tot_final_area_px_sim, tot_final_nb_cells_sim, tot_final_area_px_obs, tot_final_nb_cells_obs] = deal(zeros(7, nb_pos));
for zz = 1:nb_pos

    Position = append('/', positions(zz).name, '/');
    
    disp(Position)
    
    T_max = 10;%8;%20;%24; %Final time step in hours
    Name_file = strcat(init_name, 'fact_Pos_', num2str(zz));

    %Name used to save the dataFile (include the path, it is suggested to use the folder "Data" ):     
    DataFile = strcat(path_init,'Data/', Name_file, '.mat');
    Data_nb_cells = readtable(strcat(path_Data_position, Position, 'manual_reclass.csv'));
    
    fact_con_Pix_um = 0.0649; %Conversion factor from pixels to um. 1pixel = 0.0649um
    % mu_max_data = readtable(strcat('/Volumes/RECHERCHE/FAC/FBM/DMF/jvanderm/sinergia/D2c/UserData/Isaline/Model/', 'median_species_growth_rates_minisyncom_SE.csv')); 
    %if taking the actual measured growth rates in the SynComs
    
    mu_max_data = readtable(strcat('/Volumes/RECHERCHE/FAC/FBM/DMF/jvanderm/sinergia/D2c/UserData/Isaline/Model/', 'mean-monoculture-mu.csv')); %if taking the monoculture growth rates
    
    lag_time = load(strcat('/Volumes/RECHERCHE/FAC/FBM/DMF/jvanderm/sinergia/D2c/UserData/Isaline/Model/', 'lag_species')); 
    lag_time = lag_time.lag_species;
    lag_time = Mean_Cell(lag_time); %We will pass the lag time as Data for the SDEsSpat3D.m function
    area_species = readtable(strcat('/Volumes/RECHERCHE/FAC/FBM/DMF/jvanderm/sinergia/D2c/UserData/Isaline/Model/', 'mono_culture_kinetic_data_copy.csv')); 
    area_species.Var1 = extractBefore(area_species.Var1, 4);
    name_species_tot = unique(area_species.Var1);
    yield = load(strcat('/Volumes/RECHERCHE/FAC/FBM/DMF/jvanderm/sinergia/D2c/UserData/Isaline/Model/', 'tot_area_species')); 
    yield = yield.tot_area_species;
    yield = Mean_Cell(yield);
    max_yield = max(yield(:,1));
    yield(:,1) = yield(:,1)./max_yield*0.45; 
    %We assume the species with the maximum final area to correspond to a 45% of yield
    yield(:,2) = yield(:,2)./max_yield*0.45;
    % yield(3,:) = 3*yield(3,:);

    
    Data_kept_Pos = readtable(strcat(path_Data_position, Position, 'manual_reclass.csv')); 
    Data_kept_Pos = Data_kept_Pos(strcmp(Data_kept_Pos.channel, 'yes'), :);
    Data_nb_cells_kept = Data_nb_cells(strcmp(Data_nb_cells.channel, 'yes'), :);
    name_species = unique(Data_kept_Pos.species);
    name_species = name_species(strcmp(name_species, 'other') == 0); %Do not consider the 'other' species
    nb_species = length(name_species);
    Data_kept_Pos = Data_kept_Pos(strcmp(Data_kept_Pos.species, 'other')  == 0, :);
    name_ind = name_species;%unique(area_species.Var1);
    size_name = size(name_ind);
    mean_area = zeros(nb_species, 2);
    [~, id_name_intersect] = intersect(name_species_tot, name_species);
    yield = yield(id_name_intersect, :);
    for i = 1:size_name(1)
        ind_row = strcmp(area_species.Var1, name_ind(i));
        val_temp = table2array(area_species(ind_row, 77:84));
        mean_area(i, 1) = mean(val_temp(:), 'omitnan');
        mean_area(i, 2) = std(val_temp(:), 'omitnan');
    end
    nb_res = 6;%Change according to the matrices size below
    Pos = {};
    S_0 = zeros(1, nb_species);
    [nb_init_cells, nb_col] = size(Data_kept_Pos);
    Pos_S_0 = zeros(2, nb_init_cells);
    mu_max = zeros(nb_species, 1);
    j = 1;
    final_area_px_obs = zeros(1, nb_species);
    final_nb_cells_obs = zeros(1, nb_species);
    %Modifications here
    vect_name = []; yields_vect = [];
    for i = 1:(nb_species)
        data_temp = Data_kept_Pos(strcmp(Data_kept_Pos.species, name_species(i)), :);
        data_temp_nb_cells = Data_nb_cells_kept(strcmp(Data_nb_cells_kept.species, name_species(i)), :);
        final_area_px_obs(i) = sum(data_temp.final_area_px); %Why doesn't correspond to other measured data?
        final_nb_cells_obs(i) = sum(data_temp_nb_cells.max_cell_nrs); %Why doesn't correspond to other measured data?
        temp_x = data_temp.position_x;
        temp_y = data_temp.position_y;
        S_0(i) = length(temp_x);
        Pos_S_0(1, j:(j + S_0(i) - 1)) = fact_con_Pix_um*temp_x';
        Pos_S_0(2, j:(j + S_0(i) - 1)) = fact_con_Pix_um*temp_y';%
        j = (j + S_0(i));
        ind_mu_max = strcmp(mu_max_data.species, name_species{i});
        mu_max(i) = mu_max_data.mean_mu(ind_mu_max);%mu_max_data.medianMu_max_h_1_(ind_mu_max);
        vect_name = [vect_name; repmat(name_species(i), S_0(i), 1)];
        yields_vect = [yields_vect; repmat(yield(i,1), S_0(i), 1)];
    end
    nb_cells = length(vect_name);
    init_Dist = distEuclid(Pos_S_0, Pos_S_0); [x_ind_inf_15um, y_ind_inf_15um] = find(init_Dist < 15^2);
    remove_diag = x_ind_inf_15um - y_ind_inf_15um; 
    ind_inf = [x_ind_inf_15um(remove_diag ~= 0), y_ind_inf_15um(remove_diag ~= 0)];
    Dist_temp = zeros(nb_cells); Weighted_Mat = zeros(nb_cells);
    Dist_temp(init_Dist < 15^2) = init_Dist(init_Dist < 15^2);
    %For each cells present in the experiment, check neighbour cells from
    %different species whose distance is smaller than a defined value. Then
    %from observed value, get the factor of modification of the mu_max. If
    %there are more than one neighbour of different species present,
    %average the factors.

    %Parameters to enter
    %Create automatic number of species
    Obstacle = 0; %Species 8. 0 no obstacle, 1 there is one obstacle. The function for diffusion computation has to be changed in SDEsSpat if there are obstacles.
    dim_Img = [0 133; 0 133]; %Dimension of the frame in um
    ratio_frame = (133/dim_Img(1, 2))*(133/dim_Img(2, 2));
    
    %Attribute cell area values according to the species
    %Mass_S_0 = normrnd(2.78*10^(-13), 1e-02*2.78*10^(-13), nb_species, 1);%Mass cell at stationary phase in g.
    height_cell = normrnd(0.8, 1e-03*0.7160, nb_species, 1);
    %Height in um, corresponds to 2d where d is the circle radius. Estimated.
    length_cell = min((mean_area(:, 1)./(pi*(height_cell/2).^2) - 4/6*height_cell), 4);
    length_cell = length_cell - height_cell; 
    V_cell = pi*(height_cell/2).^2.*length_cell + 4/3*pi*(height_cell/2).^3;
    % dry weights in g C of the cells. From Appl Environ. Microbiol 64, 688–694 (1998).   
    Mass_S_0 = (421*power(V_cell(:,1),0.86)*10^(-15))./2; %(435*power(V_cell(:,1),0.86)*10^(-15))./2;     
    %Rectangle length in um. From measured microscopy data  
    Mass_Vol = Mass_S_0./(pi*(height_cell/2).^2.*(length_cell + 4/6*height_cell)); %Mass by unit of volume 


    [fact_mod_mu, fact_mod_mu_2, nb_appears, ratio_init_biomass] = deal(zeros(nb_cells, 1));
    %Only iterate on the fisrt column because of the symmetry of the pair
    for i = 1:length(ind_inf)
        species_to_mod_1 = vect_name(ind_inf(i, 1)); species_to_mod_2 = vect_name(ind_inf(i, 2));
        M0_species_1 = Mass_S_0(strcmp(species_to_mod_1, name_species)); M0_species_2 = Mass_S_0(strcmp(species_to_mod_2, name_species));
        if strcmp(species_to_mod_1, species_to_mod_2) == 0
            comb_1 = strcat(species_to_mod_1, '_', species_to_mod_2); %comb_2 = strcat(species_to_mod_2, '_', species_to_mod_1); 
            idx_1 = strcmp(table2cell(mu_Mody(:, 1)), comb_1); %idx_2 = strcmp(table2cell(mu_Mody(:, 1)), comb_2);
            % fact_mod_mu(ind_inf(i, 1)) = fact_mod_mu(ind_inf(i, 1)) + mu_Mody{idx_1, 2};
            fact_mod_mu(ind_inf(i, 2)) = fact_mod_mu(ind_inf(i, 2)) + mu_Mody{idx_1, 2};%Target that impact focal
            %fact_mod_mu(ind_inf(i, 2)) = fact_mod_mu(ind_inf(i, 2)) + mu_Mody{idx_2, 2};
            %nb_appears(ind_inf(i, 1)) = nb_appears(ind_inf(i, 1)) + 1; 
            nb_appears(ind_inf(i, 2)) = nb_appears(ind_inf(i, 2)) + 1;
            Weighted_Mat(ind_inf(i, 2), ind_inf(i, 1)) = mu_Mody{idx_1, 2};
            ratio_init_biomass(ind_inf(i, 2)) = (M0_species_1 - M0_species_2)/(1.3275e-09);
        end
    end
    for i = 1:nb_cells
        temp = sqrt(Dist_temp(i, :));
        weigth = temp(temp ~= 0)/sum(temp);
        fact_mod_mu_2(i) = sum(Weighted_Mat(i, temp ~= 0).*weigth);
    end
    fact_mod_mu = fact_mod_mu./nb_appears;
    fact_mod_mu_2(isnan(fact_mod_mu)) = 1;%For mu_max
    fact_mod_mu = fact_mod_mu_2;
    fact_mod_mu = fact_mod_mu.*yields_vect + 0*ratio_init_biomass;
    nb_cells_stat_obs = 0; %Not valid
    % fact_mod_mu = ones(nb_cells, 1); %If no variation per cell
    
    
    Threshold_divide = [normrnd(2, 1e-04, nb_species, 1) normrnd(0.2, 1e-05, nb_species, 1)];
    
    %Threshold to cell division, threshold on length. When the difference between stationary mass and the present mass is bigger or equal to this value, cell divides 
    
    %transform it into a threshold for time
    Threshold_Res = [0 0 0 0 0 0; 0 0 2 0 0 0; 0 0 5 0 0 0; 0 0 2 0 0 0; 0 0 0 0 0 0; 0 0 2 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0]; 
    
    %Threshold for resource concentration in g/mL (either put a concentration or a time threshold). %When the concentration around a cell is below this threshold, cell stop dividing. 
    %Rows equal the total number of resources and columns equal the number of species.
    std_mu_max = zeros(nb_res, nb_species);
    std_mu_max(1, :) = unifrnd(0, 0.05, 1, nb_species);%std_mu_max([1, 3], :) = unifrnd(0, 0.05, 2, nb_species);
    %Standard deviation for the mu_max (columns correspond to the species, row to the resources).
    
    kappa_3 = mu_max./yield(:,1) - mu_max;
    %Approximation obtained by the number of cells in a frame at stationary phase in comparison to the total amount of carbon resource in a frame at t = 0. 
    
    Mu_max_Mat = zeros(nb_species, nb_res); Mu_max_Mat(:, 1) = mu_max;
    ratio_mu_max = 1;
    
    %Mu_max matrix. 
    %Column for the resources, row for the species. Only the resources with a mu_max > 0 can be consumed.
    Prod_Mat = Init_Prod_Mat(kappa_3, nb_species, nb_res, 2);
    % Prod_Mat = {[0 kappa_3(1) 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];...
    %     [0 kappa_3(2) 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];...
    %     [0 kappa_3(3) 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];...
    %     [0 kappa_3(4) 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];...
    %     [0 kappa_3(5) 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];...
    %     [0 kappa_3(6) 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];...
    %     [0 kappa_3(7) 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0]}; 
    
    % Prod_Mat = {[0 kappa_3(1) 0 0 0 0; 0 0 0 0 0 0; 0 kappa_3(1) 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];...
    %     [0 kappa_3(2) 0 0 0 0; 0 0 0 0 0 0; 0 kappa_3(2) 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];...
    %     [0 kappa_3(3) 0 0 0 0; 0 0 0 0 0 0; 0 kappa_3(3) 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];...
    %     [0 kappa_3(4) 0 0 0 0; 0 0 0 0 0 0; 0 kappa_3(4) 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];...
    %     [0 kappa_3(5) 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];...
    %     [0 kappa_3(6) 0 0 0 0; 0 0 0 0 0 0; 0 kappa_3(6) 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];...
    %     [0 kappa_3(7) 0.*kappa_3(7) 0 0 0; 0 0 0 0 0 0; 0 kappa_3(7) 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0]}; 
    
    %Production rate matrices. 
    %Each row of each matrix corresponds to a resource and the rates indicate which byproduct can be %produced from the resource.
    %Larger denominator number, more sensitive to oxygen
    
    kappa_1 = 8.4710e-07*ones(nb_species, 1);
    Uptake_Mat = Uptake_Mat_Init(mu_max, kappa_1, kappa_3, nb_species, nb_res);
    % Uptake_Mat = [(mu_max(1) + kappa_3(1))/8.4710e-07 0 0 0 0 0; %9.96e-06
    %     (mu_max(2) + kappa_3(2))/8.4710e-07 0 0 0 0 0;
    %     (mu_max(3) + kappa_3(3))/8.4710e-07 0 0 0 0 0;
    %     (mu_max(4) + kappa_3(4))/8.4710e-07 0 0 0 0 0;
    %     (mu_max(5) + kappa_3(5))/8.4710e-07 0 0 0 0 0;
    %     (mu_max(6) + kappa_3(6))/8.4710e-07 0 0 0 0 0;
    %     (mu_max(7) + kappa_3(7))/8.4710e-07 0 0 0 0 0]; 
    % Uptake_Mat = [(mu_max(1) + kappa_3(1))/8.4710e-07 0 (ratio_mu_max*mu_max(1) + kappa_3(1))/8.4710e-07 0 0 0; %9.96e-06
    %     (mu_max(2) + kappa_3(2))/8.4710e-07 0 (ratio_mu_max*mu_max(2) + kappa_3(1))/8.4710e-07 0 0 0;
    %     (mu_max(3) + kappa_3(3))/8.4710e-07 0 (ratio_mu_max*mu_max(3) + kappa_3(1))/8.4710e-07 0 0 0;
    %     (mu_max(4) + kappa_3(4))/8.4710e-07 0 (ratio_mu_max*mu_max(4) + kappa_3(1))/8.4710e-07 0 0 0;
    %     (mu_max(5) + kappa_3(5))/8.4710e-07 0 0 0 0 0;
    %     (mu_max(6) + kappa_3(6))/8.4710e-07 0 (ratio_mu_max*mu_max(6) + kappa_3(1))/8.4710e-07 0 0 0;
    %     (mu_max(7) + kappa_3(7))/8.4710e-07 0 (ratio_mu_max*mu_max(7) + kappa_3(1))/8.4710e-07 0 0 0]; 
    vect_rate = Vect_rate(Uptake_Mat, Mu_max_Mat, Prod_Mat);
    shrink_rate = zeros(1, nb_species);
    Coeff_Browian = zeros(1, nb_species);
    dx = 0.5*max(length_cell(1: end - 1)); %2*max(length_cell(1: end - 1)); %4 for larges systems %um
    dx = (dim_Img(1, 2) - dim_Img(1, 1))/(floor((dim_Img(1, 2) - dim_Img(1, 1))/dx));
    dy = (dim_Img(2, 2) - dim_Img(2, 1))/(floor((dim_Img(2, 2) - dim_Img(2, 1))/dx));
    lz = 1;%round(1000/dx); %Number of resource boxes in the z-axis
    dz = 1000/lz;
    Diff_Coeff = [5760000; 1800000; 5760000; 5760000; 5760000; 7092000]; %in (um)^2/h. Increasing this diffusion will increase the consumption of the resource %Length should equals the total number of resources + wastes. It has to correspond to the number of rows and columns of the matrix Mat_rate
    t_Diff = 1./(2*Diff_Coeff*(1/dx^2 + 1/dy^2 + 1/dz^2));%dx^2./(4*Diff_Coeff); %If Lattice Blotzmann method is used
    t_D = 2e-04; %min(t_Diff) for LB; %2e-04 for CN method; %In hours. Choose the minimum diffusion time as reference time to perform consumption and diffusion.
    Time_forces = max(floor(0.0025/t_D), 1); %Data will be saved every Time_saved*t_D.
    Time_saved = Time_forces*floor((1/3)/(Time_forces*t_D));%Nb of computational iterations between two saving
    T_fin = T_max;
    vect_species = zeros(1, sum(S_0)); %Change the initialization
    ind = 1;
    for i = 1:nb_species
        vect_species(ind: ind + S_0(i) - 1) = i*ones(1, S_0(i)); %[vect_species i*ones(1, S_0(i))];
        ind = ind + S_0(i);
    end
    vect_species = [vect_species (i + 1)*ones(1, Obstacle)];
    mat_Pred = zeros(nb_species); %Predation matrix, lines correspond to predator, columns to preys %LOOK AT THIS
    Stat_Area = zeros(length(vect_species), 1);
        
    %Inital positions of the cells and initial values for the resources 
    [X, Y] = meshgrid((dim_Img(1, 1) + dx/2):dx:(dim_Img(1, 2)), (dim_Img(2, 1) + dy/2):dy:(dim_Img(2, 2)));
    Vol_Box = (1.77*10^(-8)/ratio_frame)*10^(12)/(width(X)*height(X)*lz);%(1.77*10^(-8)/ratio_frame)*10^(12)/(length(X)*length(Y)*lz);%(dim_Img(1,2)*dim_Img(2,2)*1000)*1e-18*1e03*1e12/(length(X)*length(Y)*lz);%Number of picoliter per resource box.
    
    Mass_R_0 = 7.5*1e-05*10^(-9)*Vol_Box; 
    %Concentration of 7.5e-03g/mL of soil extract in gram per resource box.
    
    %Mass_R_0 = 12.5*1e-05*10^(-9)*Vol_Box; 
    %Concentration of 12.5e-03g/mL of soil extract in gram per resource box.
    
    %2.4*10^(-4)*10^(-9)/5*Vol_Box; 
    %Concentration of 1mM of succinate in gram per resource box.
    
    %Oxygen impact is not considered
    %Mass_Oxy = 2e-04*1e-09*(Vol_Box*lz);%2*1.6*10^(-5)*10^(-9)*(Vol_Box*lz); %9e-12; %Concentration of 9mg/mL oxygen in gram per resource box on the entire length (no division by lz).
    
    Mass_R_W_0 = [Mass_R_0; 0; 0; 0; 0; 0]; %Length should equal the total number of resources + wastes. It has to correspond to the number of rows and columns of the matrix Mat_rate
    Pos_R_0 = [reshape(X, 1,[]); reshape(Y, 1,[])]; 
    Pos_Obstacle = [50; 20];%[unifrnd(dim_Img(1, 1),dim_Img(1, 2), 1, Obstacle); unifrnd(dim_Img(2, 1), dim_Img(2, 2), 1, Obstacle)];%
    
    length_cell = [length_cell, mean_area(:, 2)]; %Addition to the std for length
    [~, ~, Pos_S, vect_Cell_Length_tot, Mass_Cell_Evol, vect_angle_tot, ~, num_col, Generation_tree, Mass_Res_Waste_Evol, rho_vect, ~, lag_time, bbRegion, ~, rho_3D_tot, Vect_Time_saved] = ...
        SDEsSpat3D(lag_time, vect_species, T_fin, t_D, Diff_Coeff, t_Diff, dx, dy, Time_forces, Time_saved, length_cell, height_cell, Mass_S_0, Mass_R_W_0, Mass_Vol, Vol_Box, lz,...
        Threshold_divide, Threshold_Res, vect_rate, std_mu_max, fact_mod_mu, shrink_rate, mat_Pred, Coeff_Browian, Pos_R_0, Pos_S_0, dim_Img);
    time_step = 0:(Time_saved*t_D):T_fin;
    k = floor((1/3)/(Time_saved*t_D)); %To consider values every 20 minutes
    lag_fin = length(time_step);
    N_max = floor(lag_fin/k);
    opts_1 = odeset('RelTol',1e-9,'AbsTol',1e-10); %To smooth the curves obtained using ode45.
     
    
    
    l = 1; %Index for Biomass comparison
    [mu_max_Log_NoLT_tot, mu_max_Log_LT_tot, mu_max_Monod_LT_tot] = deal(cell(1, max(vect_species)));
    index_species = unique(vect_species);
    nb_species = length(index_species) - (Obstacle > 0);
    biomass_cell_tot = zeros(N_max,1);
     
    figure(1)
    ind_Res_Cons = 1;
    [mean_mu_max, Biomass_tot_per_species, nb_cells_stat_sim_tot, final_area_px_sim, final_nb_cells_sim] = deal(zeros(nb_species, 1));
    h = gobjects(length(name_species),1); 
    fig1 = figure(1); clf(fig1)
    ax1 = axes(fig1); hold(ax1,'on')
    for j = 1:nb_species
        Index_species = index_species(j); %Species index
        Res_Cons_ind = 1; %Resource index in order to determine the maximum growth rate associated to this resource. We are only interested in the variation of growth rates for succinate (resource 1).
        nb_col = length(num_col{j,1});
        [Lag_time_cells, mu_max_exp, mu_max_Log_NoLT, mu_max_Log_LT, mu_max_Monod_LT, Generation_Time, Doubling_Time] = deal(zeros(nb_col,1));
        scatter_Mu_max = reshape(repmat([1 2 3], nb_col, 1), [],1);
        Dist_Edge = Pos_S{j,1};
        Biomass_Edge = zeros(2, nb_col); %To asses the impact of the distance to the oxygen source on the stationary micro-colony size.
        %Iteration on the number of micro-colonies of each species
        Biomass_temp_species = 0;
        nb_cells_stat_sim_temp = 0;
        for m = 1:nb_col
            [Evol_mass,S_0, Time_step, Evol_mass_R] = deal(zeros(N_max,1));
            vect_rate_Res = vect_rate{Index_species};
            vect_rate_Res = vect_rate_Res(Res_Cons_ind, :);
            for i = 1:N_max
                temp = Mass_Cell_Evol{j, (i-1)*k + 1};
                temp_R = Mass_Res_Waste_Evol{ind_Res_Cons(1), (i-1)*k + 1}; %Change the resource index if other resource or multiple consumption. %change it according to the resource consumed
                index_col = (num_col{j, (i-1)*k + 1} == m);
                S_0(i) = sum(index_col);
                Evol_mass(i) = sum(temp(index_col));
                Evol_mass_R(i) = sum(temp_R);%Total absolute biomass in the frame
                Time_step(i) = time_step((i-1)*k + 1);
                biomass_cell_tot(i) = biomass_cell_tot(i) + Evol_mass(i);
            end
            Biomass_temp_species = Biomass_temp_species + sum(Evol_mass(end));
            Biomass_Edge(:,m) = [Dist_Edge(2,m); Evol_mass(end)];
            Lag_time_cells(m) = find(S_0 == 1, 1, 'last') - 1;

            final_nb_cells_sim(j) = final_nb_cells_sim(j) + S_0(end);
        
            Biomass_px_MeanDist = S_0*(pi*(height_cell(Index_species)/2)^2 + height_cell(Index_species)*length_cell(Index_species, 1))/(fact_con_Pix_um^2);%(Evol_mass/Mass_Vol(Index_species))/(fact_con_Pix_um^2); %Area in pixels
            final_area_px_sim(j) = final_area_px_sim(j) + Biomass_px_MeanDist(end);
            temp_LT = round(lag_time(m)/(k*Time_saved*t_D)) + 1;
            p = polyfit(Time_step, log(Biomass_px_MeanDist), 1);
            x0 = [1, 1e15*2*max(Evol_mass), 1e15*Evol_mass(1)];%[mu_max, K, R_0];
            fun_Logistic = @(x, Time_step) x(2)./(1 + (x(2) - x(3))/x(3).*exp(-x(1)*Time_step));%x(2)./(1 + exp(-x(1)*(Time_step - x(3))));%x(2)./(1 + exp(4*x(1)/x(2)*(x(3) - Time_step) + 2));%
        
            Evol_mass = 1e15*Evol_mass; %Conversion into fg
            Evol_mass_R = 1e15*Evol_mass_R; %Conversion into fg
            lb = [0, 0, 0];
            ub = [2, 20*max(Evol_mass), 20*max(Evol_mass)];
        
            r_1 = lsqcurvefit(fun_Logistic, x0, Time_step, Evol_mass, lb, ub) ; %To test with explicit solution
        
            x0 = [0.5, 2*max(Evol_mass)];
            [pars, ~] = fminsearch(@(x) -logLikelihood(Evol_mass, 1e02, Time_step, 0, x(1), lag_time(m), x(2), opts_1),...
                x0,...%Initial guess
                optimset('MaxFunEvals', 20000));
        
            sol = ode45(@(t,x) dfun(t, x, Evol_mass(1), 0, pars(1), lag_time(m), pars(2)), [0 max(Time_step)], Evol_mass(1), opts_1);
            x_est = deval(sol, Time_step);
          
            mu_max_exp(m) = p(1); mu_max_Log_NoLT(m) = r_1(1); mu_max_Log_LT(m) = pars(1); %mu_max_Monod_LT(m) = pars_2(1);
    
    
            Generation_Time(m) = log(Evol_mass(end)/Evol_mass(1) - Evol_mass(1)/Evol_mass(1))/mu_max_Monod_LT(m);
            Doubling_Time(m) = log(2)/mu_max_Monod_LT(m);
            Stat_Area(l) = Evol_mass(end); %Stationary micro-colony biomass in the order of vect_species
            l = l + 1;
            h(j) = plot(ax1, time_step(1:N_max), Evol_mass, '-', 'Color', color_set(j,:));
            hold on
	        Excel_data(:,l) = Evol_mass;
        end
        Biomass_tot_per_species(j) = Biomass_temp_species;
        mu_max_Log_NoLT_tot{j} = mu_max_Log_NoLT; mu_max_Log_LT_tot{j} = mu_max_Log_LT; mu_max_Monod_LT_tot{j} = mu_max_Monod_LT; %Growth rate for the different models
        mean_mu_max(j) = mean(mu_max_Log_LT);
        [val_sorted_Biomass_Edge, ind_sorted_Biomass_Edge] = sort(Biomass_Edge(1,:));
        scatter(ax2, final_area_px_obs(j), final_area_px_sim(j), 20, color_set(j,:), 'o', 'filled');
        scatter(ax3, final_nb_cells_obs(j), final_nb_cells_sim(j), 20, color_set(j,:), 'o', 'filled');
    end
    ratio_area_px = final_area_px_obs'./final_area_px_sim; %Ratio obs/sim
    tot_final_area_px_sim(:, zz) = final_area_px_sim; tot_final_area_px_obs(:, zz) = final_area_px_obs;
    tot_final_nb_cells_sim(:, zz) = final_nb_cells_sim; tot_final_nb_cells_obs(:, zz) = final_nb_cells_obs;
    disp(ratio_area_px)
    legend(ax1, h, name_species, 'Orientation','vertical', 'Location','southeast')
    title(ax1, 'Biomass variation with distance')
    title_str = strcat('Biomass variation with distance');
    title(title_str)
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    FigHandle = FigList(1);
    FigName = strcat(path_init, 'Figures/Fig', Name_file, num2str(1)); 
    FigName = strcat(FigName, Name_file);
    set(0, 'CurrentFigure', FigHandle);
    saveas(FigHandle, FigName, 'pdf');
    
    disp(sum(final_nb_cells_sim))
    disp(sum(final_nb_cells_obs)) 
    
    save(DataFile,'mean_mu_max', 'mu_max_Log_LT_tot', 'final_area_px_obs', 'final_area_px_sim', 'final_nb_cells_obs', 'final_nb_cells_sim', 'length_cell', 'height_cell', 'dx', 'dy', 'dim_Img',...
        'T_fin', 't_D', 'Time_saved', 'Time_forces', 'vect_species', 'Pos_S', 'vect_Cell_Length_tot', 'vect_angle_tot', 'num_col', 'Mass_Res_Waste_Evol', 'rho_vect', 'lag_time', 'bbRegion',...
        'rho_3D_tot', 'Vect_Time_saved', 'Generation_tree', 'Stat_Area')
    
    filename = append(Name_file, '_Excel','.csv');
    
    cd('/Volumes/RECHERCHE/FAC/FBM/DMF/jvanderm/sinergia/D2c/UserData/Isaline/Model/Spatial scripts/Data')
    
    writematrix(Excel_data, filename);
end
Comp_Sim_Obs_Stat_biomass = [tot_final_area_px_sim, tot_final_area_px_obs,  tot_final_nb_cells_sim, tot_final_nb_cells_obs];
T = array2table(Comp_Sim_Obs_Stat_biomass, ...
    'VariableNames', {'Area_Sim_1', 'Area_Sim_2','Area_Sim_3','Area_Sim_4', 'Area_Sim_5', 'Area_Sim_6', 'Area_Sim_7', 'Area_Sim_8',...
    'Area_Obs_1', 'Area_Obs_2','Area_Obs_3','Area_Obs_4', 'Area_Obs_5', 'Area_Obs_6', 'Area_Obs_7', 'Area_Obs_8',...
    'Cells_Sim_1', 'Cells_Sim_2','Cells_Sim_3','Cells_Sim_4', 'Cells_Sim_5', 'Cells_Sim_6', 'Cells_Sim_7', 'Cells_Sim_8',...
    'Cells_Obs_1', 'Cells_Obs_2','Cells_Obs_3','Cells_Obs_4', 'Cells_Obs_5', 'Cells_Obs_6', 'Cells_Obs_7', 'Cells_Obs_8'});
writetable(T, strcat(Data_path, init_name, 'Comp_Sim_Obs_Stat_biomass.xlsx'));