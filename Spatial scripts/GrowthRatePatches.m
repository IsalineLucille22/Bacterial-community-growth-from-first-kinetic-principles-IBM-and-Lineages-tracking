clear;
close all;


%%%%%Initialization cells positions based on real data
addpath('/Users/iguex/Library/CloudStorage/OneDrive-UniversitédeLausanne/Images_analysis/Scripts/Lineage_tracking')
Data_kept = readtable(strcat('../Data/','final_merged_table_pos1.xlsx'), 'Sheet', 1,'Format','auto'); 
Time_step = unique(Data_kept.Timepoint)/3;
Data_kept.orientation = pi/2 - Data_kept.orientation;
Data_kept.Centroid_y = 2050 - Data_kept.Centroid_y;
mean_d =  0.0649*mean(Data_kept.axis_minor_length(~isnan(Data_kept.axis_minor_length)));
mean_length =  mean(Data_kept.axis_major_length(~isnan(Data_kept.axis_major_length))) - mean_d;
Data_kept.area = Data_kept.area*0.0649^2;
Mass_S_0 = 2.78*10^(-13);
Mass_Vol = Mass_S_0./(pi*(mean_d/2)^2 + mean_d*mean_length); %Mass by unit of area
[nb_cells, mass_cells, mass_cells_2] = deal(zeros(length(Time_step) - 2, 1));
for i = 1:(length(Time_step) - 2)
    ind_t_1 = find(Data_kept.Timepoint == i);
    nb_cells(i) = length(ind_t_1);
    mass_cells_2(i) = nb_cells(i)*2.78*10^(-13); %Temporally 
    Data_1 = Data_kept(ind_t_1, :);
    ind_t_2 = find(Data_kept.Timepoint == i + 1);
    Data_2 = Data_kept(ind_t_2, :);
    Pos_X_1 = [Data_1.Centroid_x, Data_1.Centroid_y]'; %Conversion pixels into um. 1px = 0.0649 um
    [Data_kept, Pos_X_1] = Fill_NaN_Pos(ind_t_1, Data_kept, Pos_X_1);%Only for Bouke's data because of the different structure
    Pos_X_2 = [Data_2.Centroid_x, Data_2.Centroid_y]'; %Conversion pixels into um. 1px = 0.0649 um
    [Data_kept, Pos_X_2] = Fill_NaN_Pos(ind_t_2, Data_kept, Pos_X_2);%Only for Bouke's data because of the different structure
    Dist = distEuclid(Pos_X_1, Pos_X_2);
    vect_Cell_length_1 = Data_1.axis_major_length; %Conversion pixels into um. 1px = 0.0649 um
    height_cell_1 = Data_1.axis_minor_length; %Conversion pixels into um. 1px = 0.0649 um
    cell_angle_1 = pi/2 - Data_1.orientation;
    vect_Cell_length_2 = Data_2.axis_major_length; %Conversion pixels into um. 1px = 0.0649 um
    height_cell_2 = Data_2.axis_minor_length; %Conversion pixels into um. 1px = 0.0649 um
    cell_angle_2 = pi/2 - Data_2.orientation; %Conversion pixels into um. 1px = 0.0649 um
    area_opp = Data_2.area;
    %Fill the NaN values but for area
    [Data_kept, vect_Cell_length_1, height_cell_1, cell_angle_1] = Fill_NaN(ind_t_1, Dist', Data_kept, vect_Cell_length_1, height_cell_1, cell_angle_1, Pos_X_1, vect_Cell_length_2, height_cell_2, cell_angle_2, area_opp);
    mass_cells(i) = sum(Data_kept.area(ind_t_1)*Mass_Vol);
end
Evol_mass = 1e15*mass_cells; %Conversion into fg
Evol_mass_R = 5*flip(Evol_mass); %Conversion into fg
lag_time = 0;
Ks = 8.4710e+08;
vol_frame = 17700;
opts_1 = odeset('RelTol',1e-9,'AbsTol',1e-10); %To smooth the curves obtained using ode45.

%%%%%%%%%%%%
%Fit for exponential growth
%%%%%%%%%%%%

p = polyfit(Time_step(1:25), log(Evol_mass(1:25)), 1);
f_exp = @(x, t) Evol_mass(1)*exp(x*t);

f_exp_eval = f_exp(p(2), Time_step(1:20));

x0 = [0.5, 0.1, 4];%[0.5, 0.1, kappa_1, 1];%, Ks];%
[pars_2, ~] = fminsearch(@(x) -logLikelihood_2(Evol_mass, 1e01, Time_step(1: end - 2), Evol_mass_R, x(1), x(2), lag_time, Ks, x(3), 5, vol_frame, opts_1),...
x0,...%Initial guess
optimset('MaxFunEvals', 20000));

sol = ode45(@(t,x) dfun_2(t, x, pars_2(2), pars_2(1), lag_time, Ks, vol_frame, pars_2(3), 5), [0 max(Time_step)], [Evol_mass(1) Evol_mass_R(1)], opts_1);
x_est_2 = deval(sol, Time_step);
x_est_2 = x_est_2(1,:);

figure(1)
mu_max_Monod_LT = pars_2(1);
plot(Time_step(1: end - 2), Evol_mass, 'o')
hold on 
plot(Time_step, x_est_2, 'b-')
figure(2)
plot(Time_step(1:20), log(Evol_mass(1:20)))
hold on
fplot(@(t) p(1)*t + p(2), [0 Time_step(20)], 'r-')

R_init_conc = Evol_mass_R(1)/vol_frame*1e09;
mu_max_Monod_LT*R_init_conc/(R_init_conc + Ks)


Generation_Time = log(Evol_mass(end)/Evol_mass(1) - Evol_mass(1)/Evol_mass(1))/mu_max_Monod_LT;
Doubling_Time = log(2)/mu_max_Monod_LT;
Stat_Area = Evol_mass(end); %Stationary micro-colony areas in the order of vect_species