%Simulations subsytems with fitted interspecific interactions
clear
close all

%Save or Not
save_data = 0; %1 if save, 0 otherwise

%Loadind data
%Anthony's data
% Data_kept = readtable(strcat('Data/','final_merged_table.xlsx'), 'Sheet', 1,'Format','auto');
% Data_kept = readtable(strcat('Data/','final_merged_table_pos10.xlsx'), 'Sheet', 1,'Format','auto');
%Bouke's data
Data_kept = readtable(strcat('Data/','20230329_Burkholderia_PTYG_coculture_df_all_features.xlsx'), 'Sheet', 1,'Format','auto'); 
Data_kept = Data_kept(strcmp(Data_kept.Position, 'Pos10'), :);
real_sim = 0; %1 if real, 0 if sim
%Simulated data 
num_sheet = 6;%'20230329_Burkholderia_PTYG_coculture_Pos_10';%'BoukeDistBurkPos10';%'BoukePos10Dist';%3;%'Pos14';%3;
Data_kept = readtable(strcat('Data/','Simulated_Data_v2.xlsx'), 'Sheet', num_sheet,'Format','auto'); 
[~, sorted_time] = sort(Data_kept.Timepoint);
Data_kept = Data_kept(sorted_time, :);
% Data_kept.orientation = -(pi/2 - Data_kept.orientation);
% Data_kept.Centroid_y = 2050 - Data_kept.Centroid_y;
time_point = unique(Data_kept.Timepoint);
Stat_t = length(time_point); %36%Stop after this stat. time
weights = [0, 0, 0];
% len_time_point = length(time_point); %Comment or not depending on the
% final wanted time
n_1 = sum(Data_kept.Timepoint == min(time_point));
Lineage_2 = 1:n_1;
Tab_fin  = [];
To_add_next = {};
max_ID = max(Data_kept.Mask_nb) + 1;
max_Lineage = max(Lineage_2);
len_time_point = Stat_t;
threshold_New_Lin = 100;
for i = 1:(len_time_point - 2)
    t_1 = time_point(i);
    t_2 = time_point(i + 1);
    t_3 = time_point(i + 2);

    ind_t_1 = find(Data_kept.Timepoint == t_1);
    Data_1 = Data_kept(ind_t_1, :);
    ind_t_2 = find(Data_kept.Timepoint == t_2);
    Data_2 = Data_kept(ind_t_2, :);
    ind_t_3 = find(Data_kept.Timepoint == t_3);
    Data_3 = Data_kept(ind_t_3, :);
    Pos_X_1 = [Data_1.Centroid_x, Data_1.Centroid_y]';
    [Data_kept, Pos_X_1] = Fill_NaN_Pos(ind_t_1, Data_kept, Pos_X_1);%Only for Bouke's data because of the different structure
    Pos_X_2 = [Data_2.Centroid_x, Data_2.Centroid_y]';
    [Data_kept, Pos_X_2] = Fill_NaN_Pos(ind_t_2, Data_kept, Pos_X_2);%Only for Bouke's data because of the different structure
    Pos_X_3 = [Data_3.Centroid_x, Data_3.Centroid_y]';
    [Data_kept, Pos_X_3] = Fill_NaN_Pos(ind_t_3, Data_kept, Pos_X_3);%Only for Bouke's data because of the different structure
    Dist = distEuclid(Pos_X_1, Pos_X_2);
    Dist_23 = distEuclid(Pos_X_2, Pos_X_3);
    Dist_13 = distEuclid(Pos_X_1, Pos_X_3);
    Dist_23 = [Dist' Dist_23];


    vect_Cell_length_1 = Data_1.axis_major_length; 
    height_cell_1 = Data_1.axis_minor_length;
    cell_angle_1 = Data_1.orientation;
    vect_Cell_length_2 = Data_2.axis_major_length; 
    height_cell_2 = Data_2.axis_minor_length;
    cell_angle_2 = Data_2.orientation;
    vect_Cell_length_3 = Data_3.axis_major_length; 
    height_cell_3 = Data_3.axis_minor_length;
    cell_angle_3 = Data_3.orientation;

    %Fill the NaN values
    %Combination t-1 and t+1 for t_2
    [Data_kept, vect_Cell_length_1, height_cell_1, cell_angle_1] = Fill_NaN(ind_t_1, Dist', Data_kept, vect_Cell_length_1, height_cell_1, cell_angle_1, Pos_X_1, vect_Cell_length_2, height_cell_2, cell_angle_2);
    [Data_kept, vect_Cell_length_2, height_cell_2, cell_angle_2] = Fill_NaN(ind_t_2, Dist_23', Data_kept, vect_Cell_length_2, height_cell_2, cell_angle_2, Pos_X_2, [vect_Cell_length_1; vect_Cell_length_3], [height_cell_1; height_cell_3], [cell_angle_1; cell_angle_3]);

    % [Data_kept, vect_Cell_length_1, height_cell_1, cell_angle_1] = Fill_NaN(ind_t_1, Dist', Data_kept, vect_Cell_length_1, height_cell_1, cell_angle_1, Pos_X_1, vect_Cell_length_2, height_cell_2, cell_angle_2);
    n_1 = height(Data_1(:, 1));
    Seg_tot_1 = arrayfun(@(x) Rect2Seg([Pos_X_1(1,x) Pos_X_1(2,x) vect_Cell_length_1(x) height_cell_1(x)], cell_angle_1(x)),1:n_1,'UniformOutput',false); %Find an alternative
    cell_ID_1 = Data_1.Mask_nb;
    
    % [Data_kept, vect_Cell_length_2, height_cell_2, cell_angle_2] = Fill_NaN(ind_t_2, Dist, Data_kept, vect_Cell_length_2, height_cell_2, cell_angle_2, Pos_X_2, vect_Cell_length_1, height_cell_1, cell_angle_1);
    n_2 = height(Data_2(:, 1)); 
    Seg_tot_2 = arrayfun(@(x) Rect2Seg([Pos_X_2(1,x) Pos_X_2(2,x) vect_Cell_length_2(x) height_cell_2(x)], cell_angle_2(x)),1:n_2,'UniformOutput',false); %Find an alternative
    cell_ID_2 = Data_2.Mask_nb;

    [Data_kept, vect_Cell_length_3, height_cell_3, cell_angle_3] = Fill_NaN(ind_t_3, Dist_23, Data_kept, vect_Cell_length_3, height_cell_3, cell_angle_3, Pos_X_3, vect_Cell_length_2, height_cell_2, cell_angle_2);
    n_3 = height(Data_3(:, 1));
    Seg_tot_3 = arrayfun(@(x) Rect2Seg([Pos_X_3(1,x) Pos_X_3(2,x) vect_Cell_length_3(x) height_cell_3(x)], cell_angle_3(x)),1:n_3,'UniformOutput',false); %Find an alternative
    cell_ID_3 = Data_3.Mask_nb;

    Lineage_1 = Lineage_2;
    Daughter_temp = zeros(n_1, 2);

    [val, Lineage_ind] = min(Dist, [], 1);
    nb_zeros = sum(val > threshold_New_Lin); %Put a threshold value, if futher away than this value, attribute a new lineage.
    %%%%
    %%%% Comment if you don't want to add supplementary mother(s) from time
    %%%% t to time t + 1
    %Could be reattributed after!! Change the position of this loop
%     unique_mother_ID = unique(Lineage_ind);
%     No_Used_Mother_ind = find(~ismember(1:length(cell_ID_1), unique_mother_ID));
%     sum_Dist_13_row = [];
%     if ~isempty(No_Used_Mother_ind)
%         [~, ID_min_Dist_13] = min(Dist_13);
%         sum_Dist_13_row = No_Used_Mother_ind(ismember(No_Used_Mother_ind, ID_min_Dist_13));
%         if  ~isempty(sum_Dist_13_row)
%             New_ID = max(Data_kept.Mask_nb) + 1:max(Data_kept.Mask_nb) + length(sum_Dist_13_row);
%             nb_to_add = length(New_ID);
%             ind_to_add = min(ind_t_1) + sum_Dist_13_row - 1;
%             Line_to_add = Data_kept(ind_to_add, :);
%             Line_to_add.Timepoint = t_2*ones(nb_to_add, 1);
%             Line_to_add.Mask_nb = New_ID'; %New ID
%             Data_kept = [Data_kept; Line_to_add]; % Addition new lines (order will be modified at the end of the main for loop)
%             Daughter_temp(sum_Dist_13_row, :) = New_ID'.*ones(nb_to_add, 2); %Attribute their own copies as daugthers to make the bond with time t + 2
%         end
%     end
    %%%%

    %%%%
    Lineage_2 = Lineage_1(Lineage_ind);
    Lineage_2(val > threshold_New_Lin) = max_Lineage + 1:max_Lineage + nb_zeros;
    Lineage_ind(val > threshold_New_Lin) = 0;
    ind_empty = [];
    for j = 1:n_1
        ind_temp = find(Lineage_ind == j);
        ID_Daughters = cell_ID_2(ind_temp);
        if length(ID_Daughters) > 2
            [~, b] = mink(Dist(j, ind_temp), 2);
            ID_Daughters = ID_Daughters(b);
        end
        if ~isempty(ID_Daughters)
            Daughter_temp(j, :) = ID_Daughters;
        else 
            ind_empty = [ind_empty, j]; %For empty mother look at t + 1 and add them if they have daugthers 
        end
    end

    %%%% To add or remove. Impose unique lineage value for each tree
%     unique_vals = unique(Lineage_ind(Lineage_ind ~= 0));
%     if ~isempty(unique_vals)
%         counts = histcounts(Lineage_ind, [unique_vals, max(unique_vals) + 1]);
%         ind_empty = [];
%         while max(counts) > 2
%             for j = 1:n_1 
%                 ind_temp = find(Lineage_ind == j);
%                 ID_Daughters = cell_ID_2(ind_temp);
%                 if length(ID_Daughters) > 2
%                     [~, b] = mink(Dist(j, ind_temp), 2);
%                     Ind_Daughters_NK = ind_temp(setdiff(1:length(ID_Daughters), b));
%                     Dist(j, Ind_Daughters_NK) = inf;
%                     ID_Daughters = ID_Daughters(b);
%                 end
%                 if ~isempty(ID_Daughters)
%                     Daughter_temp(j, :) = ID_Daughters;
%                 else 
%                     ind_empty = [ind_empty, j];
%                 end
%             end
%             [val, Lineage_ind] = min(Dist, [], 1);
%             Lineage_2 = Lineage_1(Lineage_ind);
%             nb_zeros = sum(val > threshold_New_Lin);
%             Lineage_2(val > threshold_New_Lin) = max_Lineage + 1:max_Lineage + nb_zeros;
%             Lineage_ind(val > threshold_New_Lin) = 0;
%             ind_empty = [];
%             unique_vals = unique(Lineage_ind(Lineage_ind ~= 0));
%             counts = histcounts(Lineage_ind, [unique_vals, max(unique_vals) + 1]);
%         end
%         for j = 1:n_1 
%             ind_temp = find(Lineage_ind == j);
%             ID_Daughters = cell_ID_2(ind_temp);
%             if length(ID_Daughters) > 2
%                 [~, b] = mink(Dist(j, ind_temp), 2);
%                 Ind_Daughters_NK = ind_temp(setdiff(1:length(ID_Daughters), b));
%                 Dist(j, Ind_Daughters_NK) = inf;
%                 ID_Daughters = ID_Daughters(b);
%             end
%             if ~isempty(ID_Daughters)
%                 Daughter_temp(j, :) = ID_Daughters;
%             else 
%                 ind_empty = [ind_empty, j];
%             end
%         end
%     end
    %%%%

    Daughter_temp(Daughter_temp(:, 1) - Daughter_temp(:, 2) == 0, 2) = 0; %hlp matrix
    Tab_1 = [cell_ID_1, Daughter_temp, Lineage_1', Data_1.Timepoint,  Pos_X_1', cell_angle_1, height_cell_1, vect_Cell_length_1]; %Create here the all_tracks table 
    Tab_fin = [Tab_fin; Tab_1];
    max_Lineage = max_Lineage + nb_zeros;
    [~, sorted_time] = sort(Data_kept.Timepoint);
    Data_kept = Data_kept(sorted_time, :);
%     Lineage_2 = [Lineage_2 Lineage_1(sum_Dist_13_row)];
end
Data_kept = Data_kept(Data_kept.Timepoint < Stat_t, :);
Table_fin = array2table(Tab_fin, 'VariableNames', {'Mother_ID', 'Daughter_ID_1', 'Daughter_ID_2', 'Lineage', 'Timepoint', 'Centroid_x', 'Centroid_y', 'orientation', 'axis_minor_length', 'axis_major_length'});

cost = weight_cost(Data_kept, Table_fin, weights);
%% Creation of a lineage tree 
close all

t_0 = time_point(1);
Lineage = Table_fin.Lineage;
Lineage = unique(Lineage);
Len_Lineage = length(Lineage);

for i = 1:Len_Lineage
    Data_temp = Table_fin(Table_fin.Lineage == Lineage(i), :);
    Temp_1 = [Data_temp.Mother_ID Data_temp.Daughter_ID_1];
    Temp_1(Temp_1(:, 2) == 0, :) = [];
    Temp_2 = [Data_temp.Mother_ID Data_temp.Daughter_ID_2];
    Temp_2(Temp_2(:, 2) == 0, :) = [];
    Temp = [Temp_1; Temp_2];
    [sorted_data, ind_sorted] = sort(Temp(:,1));
    Temp = Temp(ind_sorted, :);
    T = table(Temp(:, 2), ... % CellID
           Temp(:, 1), ... % ParentID (0 indicates root/initial cell)
          'VariableNames', {'CellID', 'ParentID'});

    % Remove rows where ParentID is NaN (founder cells), or keep them for a root node
    validRows = ~isnan(T.ParentID);
    parents = T.ParentID(validRows);
    children = T.CellID(validRows);
    
    % Map CellID to unique indices
    uniqueCells = unique([T.CellID; T.ParentID(~isnan(T.ParentID))]); % All unique CellIDs
    [~, parentIdx] = ismember(parents, uniqueCells);   % Map ParentID to indices
    [~, childIdx] = ismember(children, uniqueCells);   % Map CellID to indices
    
    % Create the directed graph using the mapped indices
    lineageTree = digraph(parentIdx, childIdx);
    
    % Plot the lineage tree
    figure;
    h = plot(lineageTree, 'Layout', 'layered');
    
    % Label nodes with their CellIDs instead of indices
    labelnode(h, 1:numel(uniqueCells), string(uniqueCells));
    
    title('Lineage Tree with Mapped Cell IDs');
    xlabel('Generation');
    ylabel('Cell IDs');
    iFolderName = strcat(cd, '/Figures/');
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    for iFig = 1:length(FigList)
        FigHandle = FigList(iFig);
        FigName   = num2str(get(FigHandle, 'Number'));
        FigName = strcat('Fig', FigName);
        FigName = strcat(iFolderName, FigName, 'Tree_Dist_v', num2str(num_sheet), num2str(i));
        set(0, 'CurrentFigure', FigHandle);
        saveas(FigHandle, FigName, 'pdf');
    end
    close all
end

%% Section: Find founder cells.
% Given the difficulty to establish the true founder cells at the first
% image, we look for unique mothers that can make a lineage within the
% first 5 time points

pp=1;
missing_lineage={};

%combined new_data and next_kin

%all_tracks=[new_data(:,1:4),next_kin];
all_tracks = [Table_fin.Mother_ID Table_fin.Centroid_x Table_fin.Centroid_y Table_fin.Timepoint Table_fin.Daughter_ID_1 Table_fin.Daughter_ID_2];
% x = Data_kept.Mask_nb - Table_fin.Mother_ID;
times = time_point;

%first 5 timepoints

sub = all_tracks(ismember(all_tracks(:, 4), times(1:5)),:);

missing_table = sub;

while length(missing_table) > 0.1*length(sub)
    % to collect the lineages
    % basically we try to create a path from the mother to all the descendants,
    % including those that have no further offspring
    
    missing_founder_cells=missing_table(1,1);
    time = missing_table(1,4);
    
    fcID=missing_founder_cells;
    
    %find the corresponding row and column numbers in the data file
    
    [r,c]= find(missing_table(:,1)==fcID);
    
    if ~isempty(r)
    
        cellID=fcID; %to collect the path of cellIDs that builds the lineage
        
        %restrict the loop to the row numbers that cannot surpass the length of next_kin
        
        while r<length(missing_table)
        
	        %find the corresponding next_kin ID or daughter cell IDs by the row number
        
	        dcID = missing_table(r,5:6); 
        
	        %reshape because there will be multiple IDs to search for in subsequent generations
        
	        dcID = reshape(dcID(~isnan(dcID)),1,[]); 
        
	        % remove those cellIDs that we already have in the list
        
	        dcID(ismember(dcID,cellID))=[]; 
        
	        %find the next set of motherIDs belonging to the previous set of next_kin IDs
        
	        motherID = find(sub(:,1)==dcID); 
        
	        %find the new corresponding row numbers
        
	        [r,c] = find(missing_table(:,1)==dcID); 
        
	        %go to the next time point
        
	        time = missing_table(r,4); 
        
	        %add the new IDs to the lineage path
        
	        cellID = [cellID, dcID]; 
        
	        %empty the list of daughter IDs for the next generation
        
	        % dcID = [];
        end
    
    %collect the lineage for this founder cell and start again
    end
    
    missing_lineage{pp} = cellID;
    missing_cells = setdiff(missing_table(:,1),[missing_lineage{:}]);
    missing_table = missing_table(ismember(missing_table(:,1),missing_cells),:);
    
    pp = pp+1; %add to the tracker
end

%take from every new lineage the first motherID, which forms the founder cells

founder_cells = zeros(length(missing_lineage), 1); 

for kk = 1:length(missing_lineage)
    founder_cells(kk) = missing_lineage{kk}(1);
end


%%  Section 3: Find draft lineages %%%%%%%%%%%

% we combine the updated data and the next kin in a new array

% all_tracks=[new_data(:,1:4), next_kin];
all_tracks = [Table_fin.Mother_ID Table_fin.Centroid_x Table_fin.Centroid_y Table_fin.Timepoint Table_fin.Daughter_ID_1 Table_fin.Daughter_ID_2];

% find lineages
% loop through all cells, lineage by lineage until the table has become
% smaller than 10% of its original length. Hence: 90% of cells are tracked

pp=1;
missing_lineage={};

missing_table = all_tracks;
Len_table = size(all_tracks);
while length(missing_table) > 0.1*Len_table(1) %0.1*Len_tablelength(next_kin)

missing_founder_cells=missing_table(1,1);
% times = unique(data(:,4));
time = missing_table(1,4);

% to collect the lineages
% basically we try to create a path from the mother to all the descendants,
% including those that have no further offspring

% identify the founder cellID - keep this open so the screen output tells
% where the program is

fcID = missing_founder_cells;

%add a time tracker for the list


%find the corresponding row and column numbers in the data file

[r,c]= find(missing_table(:,1)==fcID);

cellID=fcID; %to collect the path of cellIDs that builds the lineage

%restrict the loop to the row numbers that cannot surpass the length of next_kin

while r<length(missing_table) 

	%find the corresponding next_kin ID or daughter cell IDs by the row number

	dcID=missing_table(r,5:6); 

	%reshape because there will be multiple IDs to search for in subsequent generations

	dcID=reshape(dcID(~isnan(dcID)),1,[]); 

	% remove those cellIDs that we already have in the list

	dcID(ismember(dcID,cellID))=[]; 

	%find the next set of motherIDs belonging to the previous set of next_kin IDs

	motherID=find(missing_table(:,1)==dcID); 

	%find the new corresponding row numbers

	[r,c]=find(missing_table(:,1)==dcID); 

	%go to the next time point

	time = missing_table(r,4); 

	%add the new IDs to the lineage path

	cellID = [cellID, dcID]; 

	%empty the list of daughter IDs for the next generation

	dcID=[];
end

%collect the lineage for this founder cell and start again

missing_lineage{pp} = cellID;
missing_cells=setdiff(missing_table(1:length(missing_table),1),[missing_lineage{:}]);
missing_table=missing_table(ismember(missing_table(:,1),missing_cells),:);

pp=pp+1; %add to the tracker

end

%%list lineages in an array, every column is a lineage

lng = zeros(Len_table(1), length(missing_lineage)); %to collect the lineages in the array

for kk = 1:length(missing_lineage)

    %give every new lineage a number that corresponds to the next subfile in
    %the lineage cell array by using the kk-index
    
    lng(:, kk) = double(ismember(all_tracks(:,1), missing_lineage{kk}))*kk;

end

% combine all into a single lineage column where every number is a unique
% lineage

lng_tot = sum(lng, 2);


% add to the all_tracks file as the 7th column

all_tracks(:,7) = lng_tot(1:Len_table(1));


% this is the output of this part and saved as .csv file


Variables={'MotherID','GeomX','GeomY','Timepoint','DaughterID1','DaughterID2','Lineage'};
t=array2table(all_tracks,'VariableNames',Variables);
writetable(t,'tracked_cells_Matlab.csv');


%% Section: make an intermediate plot of the tracked lineages
close all

%Keep or comment according to sections that have been run
all_tracks = [Table_fin.Mother_ID Table_fin.Centroid_x Table_fin.Centroid_y Table_fin.Timepoint Table_fin.Daughter_ID_1 Table_fin.Daughter_ID_2 Table_fin.Lineage];

figH = figure;

% subplot('Position',[0.1,0.1,0.5,0.5]);

lng_nrs = unique(all_tracks(:,7)); %all lineage numbers
colors_set_lineage = distinguishable_colors(5000);
max_Time = max(Table_fin.Timepoint);

for i = 1:length(lng_nrs)

tmp = all_tracks(all_tracks(:,7) == lng_nrs(i), :);
% tmp = tmp(tmp(:, 4) == max_Time, :);

% scatter(tmp(:,2), 1800 - tmp(:,3),'.')
scatter(tmp(:,2), tmp(:,3),[], colors_set_lineage(lng_nrs(i), :),  '.')
hold on

end

title('Stationary microcolonies');
iFolderName = strcat(cd, '/Figures/');
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
    FigHandle = FigList(iFig);
    FigName   = num2str(get(FigHandle, 'Number'));
    FigName = strcat('Fig', FigName);
    FigName = strcat(iFolderName, FigName, 'Stat_Dist_v', num2str(num_sheet));
    set(0, 'CurrentFigure', FigHandle);
    saveas(FigHandle, FigName, 'pdf');
end

%% Movie generation with cell shapes

close all;

save_data = 1;

Data_file = Table_fin;%Name of the file and the name given to the video file if created
name_movie_file = num_sheet;

% Movie generation
Tot_res = Fun_MovieSDEVec_Lineage(Data_file, real_sim, save_data, name_movie_file);