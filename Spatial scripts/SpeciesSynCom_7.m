clear;
close all;


%%%%%Initialization cells positions based on real data
%%%%%
% addpath('/Users/iguex/Library/CloudStorage/OneDrive-UniversitédeLausanne/Images_analysis/Scripts/Lineage_tracking')
% addpath('/Users/isalinelucille-guex/Library/CloudStorage/OneDrive-UniversitédeLausanne/Images_analysis/Scripts/Lineage_tracking');
path_init = '/Volumes/RECHERCHE/FAC/FBM/DMF/jvanderm/sinergia/D2c/UserData/Isaline/Model/Spatial scripts/';
path_Data_position = '/Volumes/RECHERCHE/FAC/FBM/DMF/jvanderm/sinergia/D2c/RawData/Tania/Microscopy/SE_experiments/20251022_minisyncom_SE/Dimalis';
Position = '/Pos1/';
addpath(path_init)
addpath(path_Data_position)
addpath(strcat(path_init, '/Data'))
T_max = 20;%24;%15; %Final time step in hours
Name_file = 'test_40';%append('SynCom7_SE_mu_mix','-',erase(Position,'/'));
%Name used to save the dataFile (include the path, it is suggested to use the folder "Data" ):  
DataFile = strcat('Data/', Name_file, '.mat');
fact_con_Pix_um = 0.0649; %Conversions factor from pixels to um
mu_max_data = readtable(strcat('/Volumes/RECHERCHE/FAC/FBM/DMF/jvanderm/sinergia/D2c/UserData/Isaline/Model/', 'median_species_growth_rates_minisyncom_SE.csv')); 
lag_time = load(strcat('/Volumes/RECHERCHE/FAC/FBM/DMF/jvanderm/sinergia/D2c/UserData/Isaline/Model/', 'lag_species')); 
lag_time = lag_time.lag_species;
lag_time = Mean_Cell(lag_time); %We will pass the lag time as Data for the SDEsSpat3D.m function
area_species = readtable(strcat('/Volumes/RECHERCHE/FAC/FBM/DMF/jvanderm/sinergia/D2c/UserData/Isaline/Model/', 'mono_culture_kinetic_data_copy.csv')); 
area_species.Var1 = extractBefore(area_species.Var1, 4);
yield = load(strcat('/Volumes/RECHERCHE/FAC/FBM/DMF/jvanderm/sinergia/D2c/UserData/Isaline/Model/', 'tot_area_species')); 
yield = yield.tot_area_species;
yield = Mean_Cell(yield);
max_yield = max(yield(:,1));
yield(:,1) = yield(:,1)./max_yield*0.4; %We assume the species with the maximum final area to correspond to a 30% of yield
yield(:,2) = yield(:,2)./max_yield*0.4;
Data_kept_Pos = readtable(strcat(path_Data_position, Position, 'manual_reclass.csv')); 
Data_kept_Pos = Data_kept_Pos(strcmp(Data_kept_Pos.channel, 'yes'), :);
name_species = unique(Data_kept_Pos.species);
name_species = name_species(strcmp(name_species, 'other') == 0); %Do not consider the 'other' species
nb_species = length(name_species);
Data_kept_Pos = Data_kept_Pos(strcmp(Data_kept_Pos.species, 'other')  == 0, :);
name_ind = unique(area_species.Var1);
size_name = size(name_ind);
mean_area = zeros(nb_species, 2);
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
for i = 1:(nb_species)
    data_temp = Data_kept_Pos(strcmp(Data_kept_Pos.species, name_species(i)), :);
    temp_x = data_temp.position_x;
    temp_y = data_temp.position_y;
    S_0(i) = length(temp_x);
    Pos_S_0(1, j:(j + S_0(i) - 1)) = fact_con_Pix_um*temp_x';
    Pos_S_0(2, j:(j + S_0(i) - 1)) = fact_con_Pix_um*temp_y';%(2050 - temp_y)' inversion of the axis?;
    j = (j + S_0(i));
    ind_mu_max = strcmp(mu_max_data.species, name_species{i});
    mu_max(i) = mu_max_data.medianMu_max_h_1_(ind_mu_max);
end

nb_cells_stat_obs = 0; %Not valid


%Parameters to enter
%Create automatic number of species
%Remove it
Obstacle = 0; %Species 8. 0 no obstacle, 1 there is one obstacle. The function for diffusion computation has to be changed in SDEsSpat if there are obstacles.
dim_Img = [0 133; 0 133]; %Dimension of the frame in um
ratio_frame = (133/dim_Img(1, 2))*(133/dim_Img(2, 2));
%Attribute value accrding to the species
Mass_S_0 = normrnd(2.78*10^(-13), 1e-02*2.78*10^(-13), nb_species, 1);%Mass cell at stationary phase in g.
height_cell = normrnd(0.8, 1e-03*0.7160, nb_species, 1);%Height in um, corresponds to 2d where d is the circle radius.
length_cell = min((mean_area(:, 1)./(pi*(height_cell/2).^2) - 4/6*height_cell), 4);
length_cell = length_cell - height_cell; %Rectangle length in um.
Mass_Vol = Mass_S_0./(pi*(height_cell/2).^2.*(length_cell + 4/6*height_cell)); %Mass by unit of area
Threshold_divide = [normrnd(2, 1e-04, nb_species, 1) normrnd(0.2, 1e-05, nb_species, 1)];%Threshold to cell division, threshold on length. When the difference between startionnary mass and the present mass is bigger or equal to this value, cell divides 
Threshold_Res = [0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0]; %Threshold for resource concentration in g/mL (either put a concentration or a time threshold). When the concentration around a cell is below this threshold, cell stop dividing. %Rows equal the total number of resources and columns equal the number of species.
std_mu_max = unifrnd(0, 0.1, nb_res, nb_species);%Standard deviation for the mu_max (columns correspond to the species, row to the resources).
kappa_3 = mu_max./yield(:,1) - mu_max; %Approximation obtained by the number of cells in a frame at stationary phase in comparison to the total amount of carbon resource in a frame at t = 0. 0.3 correspond to the yield
Mu_max_Mat = zeros(nb_species, nb_res); Mu_max_Mat(:, 1) = mu_max;%Mu_max matrix. %Column for the resources, row for the species. Only the resources with a mu_max > 0 can be consumed.
Prod_Mat = {[0 kappa_3(1) 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];...
    [0 kappa_3(2) 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];...
    [0 kappa_3(3) 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];...
    [0 kappa_3(4) 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];...
    [0 kappa_3(5) 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];...
    [0 kappa_3(6) 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];...
    [0 kappa_3(7) 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0]}; %Production rates matrices. %Each row of each matrix corresponds to a resource and the rates indicate which byproduct can be produced from the resource.
Uptake_Mat = [(mu_max(1) + kappa_3(1))/8.4710e-07 0 0 0 0 0; %9.96e-06
    (mu_max(2) + kappa_3(2))/8.4710e-07 0 0 0 0 0;%Larger denominator number, more sensitive to oxygen
    (mu_max(3) + kappa_3(3))/8.4710e-07 0 0 0 0 0;
    (mu_max(4) + kappa_3(4))/8.4710e-07 0 0 0 0 0;
    (mu_max(5) + kappa_3(5))/8.4710e-07 0 0 0 0 0;
    (mu_max(6) + kappa_3(6))/8.4710e-07 0 0 0 0 0;
    (mu_max(7) + kappa_3(7))/8.4710e-07 0 0 0 0 0]; 
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
vect_species = []; %Change the initialization
for i = 1:nb_species
    vect_species = [vect_species i*ones(1, S_0(i))];
end
vect_species = [vect_species (i + 1)*ones(1, Obstacle)];
mat_Pred = zeros(nb_species); %Predation matrix, lines correspond to predator, columns to preys %LOOK AT THIS
Stat_Area = zeros(length(vect_species), 1);
    
%Inital positions of the cells and initial values for the resources 
[X, Y] = meshgrid((dim_Img(1, 1) + dx/2):dx:(dim_Img(1, 2)), (dim_Img(2, 1) + dy/2):dy:(dim_Img(2, 2)));
Vol_Box = (1.77*10^(-8)/ratio_frame)*10^(12)/(width(X)*height(X)*lz);%(1.77*10^(-8)/ratio_frame)*10^(12)/(length(X)*length(Y)*lz);%(dim_Img(1,2)*dim_Img(2,2)*1000)*1e-18*1e03*1e12/(length(X)*length(Y)*lz);%Number of picoliter per resource box.
Mass_R_0 = 7.5*1e-05*10^(-9)*Vol_Box; %Concentration of 7.5e-03g/mL of soil extract in gram per resource box.%2.4*10^(-4)*10^(-9)/5*Vol_Box; %Concentration of 1mM of succinate in gram per resource box.%4.2*10^(-5)*10^(-9)*(Vol_Box*lz); %Concentration of 0.5mM of 3CBA in gram per resource box on the entire length (no division by lz).
Mass_Oxy = 2e-04*1e-09*(Vol_Box*lz);%2*1.6*10^(-5)*10^(-9)*(Vol_Box*lz); %9e-12; %Concentration of 9mg/mL oxygen in gram per resource box on the entire length (no division by lz).
Mass_R_W_0 = [Mass_R_0; 0; 0; 0; 0; 0]; %Length should equals the total number of resources + wastes. It has to correspond to the number of rows and columns of the matrix Mat_rate
Pos_R_0 = [reshape(X, 1,[]); reshape(Y, 1,[])]; 
Pos_Obstacle = [50; 20];%[unifrnd(dim_Img(1, 1),dim_Img(1, 2), 1, Obstacle); unifrnd(dim_Img(2, 1), dim_Img(2, 2), 1, Obstacle)];%

length_cell = [length_cell, mean_area(:, 2)]; %Addition to the std for length
[~, ~, Pos_S, vect_Cell_Length_tot, Mass_Cell_Evol, vect_angle_tot, ~, num_col, Generation_tree, Mass_Res_Waste_Evol, rho_vect, ~, lag_time, bbRegion, ~, rho_3D_tot, Vect_Time_saved] = ...
    SDEsSpat3D(lag_time, vect_species, T_fin, t_D, Diff_Coeff, t_Diff, dx, dy, Time_forces, Time_saved, length_cell, height_cell, Mass_S_0, Mass_R_W_0, Mass_Vol, Vol_Box, lz, Threshold_divide, Threshold_Res, vect_rate, std_mu_max, shrink_rate,...
    mat_Pred, Coeff_Browian, Pos_R_0, Pos_S_0, dim_Img);
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
 
color_set = {'red', 'blue', 'yellow', 'black', 'green', 'cyan', 'magenta', 'yellow'};
figure(1)
ind_Res_Cons = 1;
Biomass_tot_per_species = zeros(nb_species, 1);
nb_cells_stat_sim_tot = zeros(nb_species, 1);
for j = 1:nb_species
    Index_species = index_species(j); %Species index
    Res_Cons_ind = 1; %Resource index in order to determine the maximum growth rate associated to this resource. We are only interested in the variation of growth rates for succinate (resource 1).
    nb_col = length(num_col{j,1});
    [Lag_time_cells, mu_max_exp, nb_cells_stat_sim, mu_max_Log_NoLT, mu_max_Log_LT, mu_max_Monod_LT, Generation_Time, Doubling_Time] = deal(zeros(nb_col,1));
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
        nb_cells_stat_sim(m) = S_0(end); %Number of cells into the microcolony at stationary phase (or last time step)
        nb_cells_stat_sim_temp = nb_cells_stat_sim_temp + nb_cells_stat_sim(m);
    
        Biomass_MeanDist = (Evol_mass/Mass_Vol(Index_species))/(fact_con_Pix_um^2);
        temp_LT = round(lag_time(m)/(k*Time_saved*t_D)) + 1;
        p = polyfit(Time_step, log(Biomass_MeanDist), 1);
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
        Stat_Area(l) = Evol_mass(end); %Stationary micro-colony areas in the order of vect_species
        l = l + 1;
        plot(1:N_max, Evol_mass, '-', 'Color', color_set{j})
        hold on
    end
    Biomass_tot_per_species(j) = Biomass_temp_species;
    nb_cells_stat_sim_tot(j) = nb_cells_stat_sim_temp;
    mu_max_Log_NoLT_tot{j} = mu_max_Log_NoLT; mu_max_Log_LT_tot{j} = mu_max_Log_LT; mu_max_Monod_LT_tot{j} = mu_max_Monod_LT; %Growth rate for the different models
    [val_sorted_Biomass_Edge, ind_sorted_Biomass_Edge] = sort(Biomass_Edge(1,:));
    Excel_data = [val_sorted_Biomass_Edge' Biomass_Edge(2,ind_sorted_Biomass_Edge)' mu_max_Monod_LT(ind_sorted_Biomass_Edge)];
    writematrix(Excel_data, strcat(path_init, 'Results/', Name_file, name_species{index_species(j)}, 'Excel_data.xlsx'))
end
legend(name_species, 'Orientation', 'vertical', 'Location', 'southeast')
title_str = strcat('Biomass variation with distance');
title(title_str)
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
FigHandle = FigList(1);
FigName = strcat(path_init, 'Figures/Fig', Name_file, num2str(1)); 
FigName = strcat(FigName, Name_file);
set(0, 'CurrentFigure', FigHandle);
saveas(FigHandle, FigName, 'pdf');

disp(sum(nb_cells_stat_sim_tot)) %Check it

save(DataFile,'length_cell', 'height_cell', 'dx', 'dy', 'dim_Img', 'T_fin', 't_D', 'Time_saved', 'Time_forces', 'vect_species', 'Pos_S', 'vect_Cell_Length_tot', 'vect_angle_tot', 'num_col', 'Mass_Res_Waste_Evol', 'rho_vect', 'lag_time', 'bbRegion',...
    'rho_3D_tot', 'Vect_Time_saved', 'Generation_tree', 'nb_cells_stat_sim', 'nb_cells_stat_obs','Stat_Area')

%% Video generation

clear;
close all;

% Add path to run the code from main folder
Position = '/Pos7/';
save_data = 1; %Create piece of frame

Name_file = 'test_40';%'SynCom7_SE_mu_mix-Pos1';%append('SynCom7_SE','-',erase(Position,'/'));%Name of the file and the name given to the video file if created
resultsFile = strcat('Data/', Name_file, '.mat'); 
res_number = 1; %Number of the different resources (resources and wastes)

% Movie generation
Delta_Time = 20;%60; %Number of minutes between two video figures. This value should be equal or higher than the saved time defined before
Tot_res = Fun_MovieSDEVec_Lineage(res_number, resultsFile, save_data, Name_file, Delta_Time);