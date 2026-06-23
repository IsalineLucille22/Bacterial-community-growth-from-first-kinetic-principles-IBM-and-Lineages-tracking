function [Tot_res] = Fun_SpecificTime_Growth_rates(Data, real_sim, save_data, specified_time_point)
time_point = unique(Data.Timepoint);
T_fin = time_point(end);
Pos_fin_x = Data.Centroid_x(Data.Timepoint == T_fin); Pos_fin_y = Data.Centroid_y(Data.Timepoint == T_fin);
min_x = min(Pos_fin_x); max_x = max(Pos_fin_x);
min_y = min(Pos_fin_y); max_y = max(Pos_fin_y);
dim_Img = real_sim*[(min_x - 100) (max_x + 100); (min_y - 100) (max_y + 100)] + (1 - real_sim)*[(min_x - 13) (max_x + 13); (min_y - 13) (max_y + 13)];
height_cell_mean = mean(Data.axis_minor_length);
length_cell_mean = mean(Data.axis_major_length) - real_sim*height_cell_mean;
x_axis = [dim_Img(1, 1) dim_Img(1, 2)];
y_axis = [dim_Img(2, 1) dim_Img(2, 2)];
p_Diff = time_point(2) - time_point(1); 
T_fin = time_point(end);
n_step = floor((T_fin - min(Data.Timepoint) + 1)/p_Diff);
mu_max_unique = Data.mu_max_vect;
%division_rate = log(2)./mu_max_unique;
%mu_max_unique = division_rate;
%limit mu_max to avoid abnormalities

mu_max_unique(mu_max_unique>1)=1;

mu_max_normalized = (mu_max_unique + 1 - min(mu_max_unique + 1)) / (max(mu_max_unique + 1) - min(mu_max_unique + 1)); % Normalize to [0, 1]
%colormap_used = jet(256);
%make shades in blue intensity

% 256-color palette: white -> blue-ish dark blue (hue held roughly constant)
n = 256;

% Choose a target “dark blue” (edit these if you want a different end color)
targetHex = '#00008B';                 % DarkBlue (0,0,139)
    hexStr = lower(targetHex);
    hexStr = hexStr(1+startsWith(hexStr,'#'):end);
    r = sscanf(hexStr(1:2), '%x');
    g = sscanf(hexStr(3:4), '%x');
    b = sscanf(hexStr(5:6), '%x');
    rgb = double([r g b]) / 255;  % [0..1] RGB from hex
targetH = rgb2hsv(rgb);       % hue/sat/value of the target

% We'll hold hue near the target's hue and interpolate saturation + lightness-ish via V.
% Approach: interpolate in HSV with hue fixed, and V decreasing from 1 to target V.
h = targetH(1);
s0 = 0;               % white has ~0 saturation
s1 = targetH(2);     % end saturation (blue-ness)
v0 = 1;               % white value/brightness
v1 = targetH(3);     % end brightness

palette = zeros(n,3);
for zz = 1:n
    t = (zz-1)/(n-1);
    si = s0 + t*(s1 - s0);
    vi = v0 + t*(v1 - v0);
    hi = h;  % keep hue fixed (blue-ish)

    % Convert HSV -> RGB
    rgb = hsv2rgb([hi, si, vi]);
    palette(zz,:) = rgb;
end

colormap_used=palette;

colors_set_lineage_col = colormap_used(max(ceil(mu_max_normalized*length(colormap_used)), 1), :);

Nb_species = 1;%length(index_species);
color_set = {'red', 'blue', 'yellow', 'black', 'green', 'cyan', 'magenta', 'yellow'};
scale_factor = 1;

fig = {};
for i = 1:Nb_species
    fig{i} = scatter(NaN, NaN, height_cell_mean^2*pi*353/scale_factor^2, color_set{i}, 'filled', 'MarkerEdgeColor', [0 0 0], 'LineWidth', 1);
    hold on
end
colormap(colormap_used); % Set the colormap for the figure
%caxis([min(mu_max_unique), max(mu_max_unique)]); % Set color axis limits based on mu_max values
caxis([0,0.6]);
colorbar; % Display the colorbar
set(gca, 'XLim', scale_factor*dim_Img(1,:), 'Ylim', scale_factor*dim_Img(2,:), 'YTickLabel', [],'XTickLabel', [], 'YTick', [], 'XTick', [], 'Visible','off')%, 'YDir','normal')
axis('image');
Tot_res = zeros(1,n_step);
rectangle2([0 0  0 0],'Curvature',[4*0.133*ones(1,1) 1*ones(1,1)], 'Rotation', 0, 'FaceColor', color_set(1), 'EdgeColor', [0 0 0], 'LineWidth', 3);
axis equal;
axis([x_axis(1) x_axis(2) y_axis(1) y_axis(2)]);

delete(findobj('type', 'patch'));
t_temp = time_point(specified_time_point);
ind_col = find(Data.Timepoint == t_temp);
Data_temp = Data(Data.Timepoint == t_temp, :);
for i = 1:Nb_species
    P_temp_S = [Data_temp.Centroid_x'; Data_temp.Centroid_y'];
    height_cell = Data_temp.axis_minor_length;
    vect_Cell_Length_temp = Data_temp.axis_major_length - real_sim*height_cell;
    vect_angle_temp = Data_temp.orientation;
    vect_angle_temp = rad2deg(vect_angle_temp);
    if ~isempty(vect_Cell_Length_temp)
       rectangle2([P_temp_S(1,:)' (P_temp_S(2,:))'  (vect_Cell_Length_temp + height_cell) height_cell],'Curvature',...
            [min((length_cell_mean + height_cell)./(vect_Cell_Length_temp + height_cell)*0.5882.*ones(length(vect_angle_temp),1), 1) 1*ones(length(vect_angle_temp),1)],...
            'Rotation', vect_angle_temp, 'FaceColor', colors_set_lineage_col(ind_col, :), 'EdgeColor', [0 0 0], 'LineWidth',0.3);
    end
end
text_time = text(0.8, 1, strcat('T :',num2str(t_temp)), 'FontSize', 6);

