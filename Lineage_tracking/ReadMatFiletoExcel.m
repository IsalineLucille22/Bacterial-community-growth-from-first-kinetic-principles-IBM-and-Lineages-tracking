clear
close all

%Save or Not
save_data = 0; %1 if save, 0 otherwise

name_Mat = 'Anthony_Lineage_v6.mat';
Excel_file = 'Data/Simulated_Data_v2.xlsx'; 
Mat_file = load(strcat('Data/', name_Mat));
Mat_file.Generation_tree{1} = Mat_file.Generation_tree{2};

Lineage_Mat_fin = fun_Lineage_Mat(Mat_file.Generation_tree);
lineage_temp = Mat_file.num_col;
lineage_temp = [lineage_temp{:}]';
Lineage_Mat_fin = [Lineage_Mat_fin, lineage_temp(1:length(Lineage_Mat_fin(:, 1)))];
Table_fin = array2table(Lineage_Mat_fin, 'VariableNames', {'Mother_ID', 'Daughter_ID_1', 'Daughter_ID_2', 'Lineage'});

Name_columns = {'Generation_tree', 'Pos_S', 'vect_Cell_Length_tot', 'vect_angle_tot', 'num_col', 'Vect_Time_saved'};
Positions = {'A2', 'B2', 'I2', 'K2', 'L2', 'D2'};
num_col = length(Name_columns);

for i = 1:num_col
    name_temp = Name_columns{i};
    data = Mat_file.(name_temp);
    data = [data{:}]';
    writematrix(data, Excel_file, 'Sheet', 6,'Range', Positions{i})
end


%% Creation of a lineage tree 
close all

Lineage = Table_fin.Lineage;
Lineage = unique(Lineage);
Len_Lineage = length(Lineage);

for i = 1:Len_Lineage
    Data_temp = Table_fin(Table_fin.Lineage == Lineage(i), :);
    [s, m] = size(Data_temp);
    if s > 1
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
            FigName = strcat(iFolderName, FigName, name_Mat, num2str(i));
            set(0, 'CurrentFigure', FigHandle);
            saveas(FigHandle, FigName, 'pdf');
        end
        close all
    end
end




