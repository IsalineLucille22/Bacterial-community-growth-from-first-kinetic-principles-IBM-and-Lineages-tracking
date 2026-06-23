clear;
close all;

Position = strcat('/Pos', string(1), '/');
path_init = '/Volumes/RECHERCHE/FAC/FBM/DMF/jvanderm/sinergia/D2c/UserData/Isaline/Model/Spatial scripts/';
path_Data_position = '/Volumes/RECHERCHE/FAC/FBM/DMF/jvanderm/sinergia/D2c/RawData/Tania/Microscopy/SE_experiments/20251022_minisyncom_SE/Dimalis';
Data_kept_Pos = readtable(strcat(path_Data_position, Position, 'manual_reclass.csv')); 
name_species = unique(Data_kept_Pos.species); name_species = name_species(1:end - 1);
nb_species = numel(name_species);
cd(path_Data_position)
positions = dir('Pos*');
nb_pos = length(positions);
cd(path_init)

fig1 = figure(1); clf(fig1)
ax1 = axes(fig1); hold(ax1,'on')
xlabel(ax1,'Simulated area')
ylabel(ax1,'Observed area')
xlim(ax1,[0 3e05])
ylim(ax1,[0 3e05])
fplot(ax1, @(x) x, [0 3e05], 'k--', 'Linewidth', 0.5)
fig2 = figure(2); clf(fig2)
ax2 = axes(fig2); hold(ax2,'on')
xlabel(ax2,'Simulated area')
ylabel(ax2,'Observed area')
xlim(ax2,[0 1e03])
ylim(ax2,[0 1e03])
fplot(ax2, @(x) x, [0 1e03], 'k--', 'Linewidth', 0.5)
color_set = [
    1 0 0;   %red
    0 0 1;   %blue
    1 0.5 0; %orange
    0 0 0;   %black
    0 1 0;   %green
    0 1 1;   %cyan
    1 0 1;   %magenta
];

h = gobjects(nb_species,1);
[val_area_per_species_obs, val_area_per_species_sim] = deal(zeros(nb_species, nb_pos));
for k = 1:nb_pos
    Pos = k;
    Name_file = strcat('SynCom7_SE_mono_yield_with_competition_index_corrections_Pos_', string(Pos));%strcat('SynCom7_SE_minisyncom_growth_rate_wo_umax_fact_Pos_', string(Pos));%'SynCom7_SE_mono_growth_rate_with_umax_fact_Pos_', 'SynCom7_SE_minisyncom_growth_rate_wo_umax_fact_Pos_'
    filename = fullfile(path_init, 'Data', Name_file + ".mat");
    load(strcat(path_init,'Data/', Name_file, '.mat'));
    for i = 1:nb_species
        h(i) = scatter(ax1, final_area_px_sim(i), final_area_px_obs(i), 40, color_set(i,:), 'o', 'filled'); %sim~obs
        scatter(ax2, final_nb_cells_sim(i), final_nb_cells_obs(i), 40, color_set(i,:), 'o', 'filled');
        val_area_per_species_obs(i, k) = final_area_px_obs(i);
        val_area_per_species_sim(i, k) = final_area_px_sim(i);
    end
end

legend(ax1, h, name_species,'Orientation','vertical', 'Location','southeast')

%%%%Linear regression, %sim~obs%%%%
[p_val, R2_LR, slope, intercept, mean_ratio] = deal(zeros(nb_species, 1));
for i = 1:nb_species
    mdl = fitlm(val_area_per_species_sim(i,:), val_area_per_species_obs(i,:));%val_area_per_species_sim(i,:) + normrnd(0, 10e-3*val_area_per_species_sim(i,:)));
    coeffs = mdl.Coefficients.Estimate;
    CovB = mdl.CoefficientCovariance;
    intercept(i) = coeffs(1);
    slope(i) = coeffs(2);
    se_slope = sqrt(CovB(2,2)); %Standard error
    t_slope_1 = (slope(i) - 1)/se_slope; %t-test h0: slope = 1
    df = mdl.DFE;
    p_val(i) = 2*(1 - tcdf(abs(t_slope_1), df));
    R2_LR(i) = mdl.Rsquared.Ordinary;
    figure(2 + i)
    plot(mdl)
    mean_ratio(i) = mean(val_area_per_species_sim(i,:)./val_area_per_species_obs(i,:));

    %%%RMA regression%%%%
    x = val_area_per_species_sim(i,:);
    y = val_area_per_species_obs(i,:);
    r = corr(x', y');

    slope_RMA = sign(r)*std(y)/std(x);
    intercept_RMA = mean(y) - slope_RMA*mean(x);
    
    y_RMA = intercept_RMA + slope_RMA*x;
    
    SS_res = sum((y - y_RMA).^2);
    SS_tot = sum((y - mean(y)).^2);
    R2 = r^2;

    
    n = nb_pos;
    
    % Approximate SE of log(slope)
    SE_log_slope = sqrt((1 - r^2)/(n - 2));
    
    alpha = 0.05;
    tcrit = tinv(1 - alpha/2, n - 2);
    
    CI_slope = slope_RMA*exp([-1 1]*tcrit*SE_log_slope);

    B = 10000;
    slope_boot = nan(B,1);
    
    for b = 1:B
        ii = randi(n, n, 1);
        xb = val_area_per_species_sim(ii);
        yb = val_area_per_species_obs(ii);
    
        rb = corr(xb, yb);
        slope_boot(b) = sign(rb) * std(yb) / std(xb);
    end
    
    CI_boot = prctile(slope_boot, [2.5 97.5]);
    
    disp(CI_boot)
end