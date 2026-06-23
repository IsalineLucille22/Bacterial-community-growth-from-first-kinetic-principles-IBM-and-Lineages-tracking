function [Tot_res] = Fun_MovieSDEVec_Lineage(Data, real_sim, save_data, name_movie_file, initial_time_point)
time_point = unique(Data.Timepoint);
T_fin = time_point(end);
Pos_fin_x = Data.Centroid_x(Data.Timepoint == T_fin); Pos_fin_y = Data.Centroid_y(Data.Timepoint == T_fin);
min_x = min(Pos_fin_x); max_x = max(Pos_fin_x);
min_y = min(Pos_fin_y); max_y = max(Pos_fin_y);
dim_Img = real_sim*[(min_x - 100) (max_x + 100); (min_y - 100) (max_y + 100)] + (1 - real_sim)*[(min_x - 13) (max_x + 13); (min_y - 13) (max_y + 13)];
height_cell_mean = mean(Data.axis_minor_length);
length_cell_mean = mean(Data.axis_major_length) - real_sim*height_cell_mean;
res_number_ind = 1;
x_axis = [dim_Img(1, 1) dim_Img(1, 2)];
y_axis = [dim_Img(2, 1) dim_Img(2, 2)];
p_Diff = time_point(2) - time_point(1); 
T_fin = time_point(end);
n_step = floor((T_fin - min(Data.Timepoint) + 1)/p_Diff);

name_movie_file = strcat(name_movie_file, num2str(res_number_ind));
Nb_species = 1;%length(index_species);
Num_color_set = {[1;3]; 1; [2;3]; 2; 3; 1; [1;2]};
color_set = {'red', 'blue', 'yellow', 'black', 'green', 'cyan', 'magenta', 'yellow'};
colors_set_lineage = distinguishable_colors(1500);
scale_factor = 1;

if save_data == 1
	myVideo = VideoWriter(strcat('Videos/',name_movie_file));
	myVideo.FrameRate = 10;
    open(myVideo)
end

fig = {};
for i = 1:Nb_species
    fig{i} = scatter(NaN, NaN, height_cell_mean^2*pi*353/scale_factor^2, color_set{i}, 'filled', 'MarkerEdgeColor', [0 0 0], 'LineWidth', 1);
    hold on
end
set(gca, 'XLim', scale_factor*dim_Img(1,:), 'Ylim', scale_factor*dim_Img(2,:), 'YTickLabel', [],'XTickLabel', [], 'YTick', [], 'XTick', [], 'Visible','off')%, 'YDir','normal')
axis('image');
t = 0;
k = (1 - initial_time_point)*(n_step - 1) + initial_time_point*1;
Tot_res = zeros(1,n_step);
rectangle2([0 0  0 0],'Curvature',[4*0.133*ones(1,1) 1*ones(1,1)], 'Rotation', 0, 'FaceColor', color_set(1), 'EdgeColor', [0 0 0], 'LineWidth', 3);
axis equal;
axis([x_axis(1) x_axis(2) y_axis(1) y_axis(2)]);

% Define unique lineages and their corresponding colors
unique_lineages = unique(Data.Lineage);
legend_colors = colors_set_lineage(unique_lineages, :); % Get colors for each unique lineage

while k < n_step
    delete(findobj('type', 'patch'));
    text_time = text(0.8*dim_Img(1,1), 0.8*dim_Img(1,1), strcat('Time (h):',num2str(t)));
    t_temp = time_point(k);
    ind_col = find(Data.Timepoint == t_temp);
    Data_temp = Data(Data.Timepoint == t_temp, :);
    for i = 1:Nb_species
        index_lineage_temp = Data_temp.Lineage;
        P_temp_S = [Data_temp.Centroid_x'; Data_temp.Centroid_y'];
        Lineage_temp = Data_temp.Lineage;
        height_cell = Data_temp.axis_minor_length;
        vect_Cell_Length_temp = Data_temp.axis_major_length - real_sim*height_cell;
        vect_angle_temp = Data_temp.orientation;
        vect_angle_temp = rad2deg(vect_angle_temp);
        if ~isempty(vect_Cell_Length_temp)
            rectangle2([P_temp_S(1,:)' (P_temp_S(2,:))'  (vect_Cell_Length_temp + height_cell) height_cell],'Curvature',...
                [min((length_cell_mean + height_cell)./(vect_Cell_Length_temp + height_cell)*0.5882.*ones(length(vect_angle_temp),1), 1) 1*ones(length(vect_angle_temp),1)],...
                'Rotation', vect_angle_temp, 'FaceColor', colors_set_lineage(index_lineage_temp, :), 'EdgeColor', [0 0 0], 'LineWidth',0.3);
               %text(P_temp_S(1, :), P_temp_S(2, :), num2str(Lineage_temp), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 6, 'Color', 'black');
        end
    end
    pause(p_Diff);
    t = time_point(k); %t + 0*p_Diff;
    k = k + 1;
    if save_data == 1
        frame = getframe(gcf);
        writeVideo(myVideo, frame);
    end
    delete(text_time);
end



if save_data == 1
    close(myVideo)

    lineage_number = unique(Data.Lineage);
    if length(lineage_number) ~= 1
        lineage_number = 'tot';
    end
    iFolderName = strcat(cd, '/Figures/');
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    for iFig = 1:length(FigList)
        FigHandle = FigList(iFig);
        FigName   = num2str(get(FigHandle, 'Number'));
        FigName = strcat('Fig', FigName);
        FigName = strcat(iFolderName, FigName, 'Stat_time_lineage', num2str(lineage_number));
        set(0, 'CurrentFigure', FigHandle);
        saveas(FigHandle, FigName, 'pdf');
    end
end