function Data_kept = fun_Slinding_Effect_v2(Data_kept)
time_point = unique(Data_kept.Timepoint);
% centroids = zeros(length(time_point), 2);
timepoints = Data_kept.Timepoint;
positions = [Data_kept.Centroid_x, Data_kept.Centroid_y];
% cell_ids = Data_kept.Mask_nb;
positions_by_time = cell(length(time_point), 1);
for i = 1:length(time_point)
    positions_by_time{i} = positions(timepoints == time_point(i), :);
end
drift = zeros(length(time_point) - 1, 2); % Drift in X and Y

for i = 1:length(time_point) - 1
    current_positions = positions_by_time{i};
    next_positions = positions_by_time{i + 1};
    pairwise_distances = pdist2(current_positions, next_positions);
    [~, closest_indices] = min(pairwise_distances, [], 2);
    matched_positions = next_positions(closest_indices, :);
    displacements = matched_positions - current_positions;
    drift(i, :) = median(displacements, 1, 'omitnan');
end
ind_mod = vecnorm(drift') <= 5; %Put the 5 as an input of the function
drift(ind_mod, :) = 0;
corrected_positions = positions;
cumulative_drift = [0, 0];
for i = 2:length(time_point)
    cumulative_drift = cumulative_drift + drift(i - 1, :);
    current_indices = (timepoints == time_point(i));
    corrected_positions(current_indices, :) = positions(current_indices, :) - cumulative_drift;
end

% Add corrected positions to the table
Data_kept.Centroid_x = corrected_positions(:, 1);
Data_kept.Centroid_y = corrected_positions(:, 2);
% for i = 1:length(time_point)
%     ind_temp = find(Data_kept.Timepoint == time_point(i));
%     centroids(i, :) = mean([Data_kept.Centroid_x(ind_temp), Data_kept.Centroid_y(ind_temp)], 1);
% end
% shifts = diff(centroids);
% correctedPositions = Data_kept; % Copy original data
% cumulativeShift = [0, 0]; % Initialize cumulative shift [X, Y]
% for i = 2:length(time_point)
%     ind_temp = correctedPositions.Timepoint == time_point(i);
%     cumulativeShift = cumulativeShift + shifts(i - 1, :);
%     correctedPositions.Centroid_x(ind_temp) = correctedPositions.Centroid_x(ind_temp) - cumulativeShift(1);
%     correctedPositions.Centroid_y(ind_temp) = correctedPositions.Centroid_y(ind_temp) - cumulativeShift(2);
% end
% Data_kept = correctedPositions;
end

% nb_cells = max(Data_kept.Mask_nb);
% nb_times = length(time_point);
% Positions_Mat = nan(nb_times, nb_cells*2);
% for i = 1:nb_times
%     temp_time = time_point(i);
%     ind_temp = Data_kept.Timepoint == time_point(i);
%     Pos_X = Data_kept.Centroid_x(ind_temp); Pos_Y = Data_kept.Centroid_y(ind_temp);
%     Positions_Mat(i, 1:2:end) = Pos_X; %Alternate the x and y positions
%     Positions_Mat(i, 2:2:end) = Pos_Y;
% end
% 
% [correctedPositions, drift] = drifty_shifty_deluxe(Positions_Mat);