clear;
close all;

Name_file_init = 'IG_SynCom7_SE_mono_yield_with_competition_index_corrections_Pos_';%'SynCom7_SE_minisyncom_growth_rate_wo_umax_fact_Pos_'; %'SynCom7_SE_mono_growth_rate_wo_umax_fact_Pos_', 'SynCom7_SE_mono_growth_rate_with_umax_fact_Pos_', 'SynCom7_SE_minisyncom_growth_rate_wo_umax_fact_Pos_'
Position = strcat('/Pos', string(1), '/');
path_init = '/Volumes/RECHERCHE/FAC/FBM/DMF/jvanderm/sinergia/D2c/UserData/Isaline/Model/Spatial scripts/';
path_Data_position = '/Volumes/RECHERCHE/FAC/FBM/DMF/jvanderm/sinergia/D2c/RawData/Tania/Microscopy/SE_experiments/20251022_minisyncom_SE/Dimalis';
path_figs_saved = strcat(path_init, 'Scripts paper IG/', 'Figures/'); 
Data_kept_Pos = readtable(strcat(path_Data_position, Position, 'manual_reclass.csv')); 
name_species = unique(Data_kept_Pos.species); name_species = name_species(1:end - 1);
nb_species = numel(name_species);
cd(path_Data_position)
positions = dir('Pos*');
nb_pos = length(positions);
cd(path_init)

fig2 = figure(2); clf(fig2)
ax2 = axes(fig2); hold(ax2,'on')
xlabel(ax2,'Simulated area')
ylabel(ax2,'Observed area')
xlim(ax2,[0 3e05])
ylim(ax2,[0 3e05])
fplot(ax2, @(x) x, [0 3e05], 'k--', 'Linewidth', 0.5)
fig3 = figure(3); clf(fig3)
ax3 = axes(fig3); hold(ax3,'on')
xlabel(ax3,'Simulated nb cells')
ylabel(ax3,'Observed nb cells')
xlim(ax3,[0 1e03])
ylim(ax3,[0 1e03])
fplot(ax3, @(x) x, [0 1e03], 'k--', 'Linewidth', 0.5)
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
    %strcat('SynCom7_SE_mono_yield_with_competition_index_corrections_Pos_', string(Pos));%
    Name_file = strcat(Name_file_init, string(Pos));%'SynCom7_SE_mono_growth_rate_wo_umax_fact_Pos_', 'SynCom7_SE_mono_growth_rate_with_umax_fact_Pos_', 'SynCom7_SE_minisyncom_growth_rate_wo_umax_fact_Pos_'
    filename = fullfile(path_init, 'Data', Name_file + ".mat");
    load(strcat(path_init,'Data/', Name_file, '.mat'));
    for i = 1:nb_species
        h(i) = scatter(ax2, final_area_px_sim(i), final_area_px_obs(i), 40, color_set(i,:), 'o', 'filled'); %sim~obs
        scatter(ax3, final_nb_cells_sim(i), final_nb_cells_obs(i), 40, color_set(i,:), 'o', 'filled');
        val_area_per_species_obs(i, k) = final_area_px_obs(i);
        val_area_per_species_sim(i, k) = final_area_px_sim(i);
    end
end

legend(ax2, h, name_species,'Orientation','vertical', 'Location','southeast')
saveas(fig3, fullfile(path_figs_saved,  strcat(Name_file_init,'simulated_vs_observed_cells.pdf')))

%%%%Linear regression, %sim~obs%%%%
[p_val_joint, R2, slope, intercept] = deal(zeros(nb_species, 1));
for i = 1:nb_species
    mdl = fitlm(val_area_per_species_sim(i,:), val_area_per_species_obs(i,:));
    coeffs = mdl.Coefficients.Estimate;
    CovB = mdl.CoefficientCovariance;
    intercept(i) = coeffs(1);
    slope(i) = coeffs(2);
    df = mdl.DFE;
    R2(i) = mdl.Rsquared.Ordinary;
    figure
    plot(mdl)

    x = val_area_per_species_sim(i,:)';
    y = val_area_per_species_obs(i,:)';
      
    H = [1 0; 0 1];    
    c = [0; 1];
    [p_val_joint(i), F] = coefTest(mdl, H, c);
    
    disp(p_val_joint(i))
end

figure(fig2)
axes(ax2)

txt = strings(nb_species,1);
for i = 1:nb_species
    txt(i) = sprintf('%s: p = %.3g, R^2 = %.2f', ...
        name_species{i}, p_val_joint(i), R2(i));
end

annotation(fig2, 'textbox', [0.15 0.55 0.35 0.35], ...
    'String', txt, ...
    'FitBoxToText', 'on', ...
    'BackgroundColor', 'white', ...
    'EdgeColor', 'black');

saveas(fig2, fullfile(path_figs_saved, strcat(Name_file_init, 'simulated_vs_observed_area.pdf')))


%%%MDS plot%%%

A = [val_area_per_species_obs val_area_per_species_sim]';%Row: replicates, columns: species

% Optional: normalize abundances to relative abundance
A_rel = A./sum(A, 2);

% Compute pairwise dissimilarity between communities
Dist = BrayCurtisDistance(A_rel);  % good for abundance/community data

% Classical MDS/PCoA
[Y, eigvals] = cmdscale(Dist, 2);

group = [repmat("Observed", nb_pos, 1); repmat("Simulated", nb_pos, 1)];

fig10 = figure;
gscatter(Y(:,1), Y(:,2), group);
xlabel('MDS1');
ylabel('MDS2');
title('MDS of final abundances');
grid on;

saveas(fig10, fullfile(path_figs_saved, strcat(Name_file_init, 'MDS_plot_Community.pdf')))