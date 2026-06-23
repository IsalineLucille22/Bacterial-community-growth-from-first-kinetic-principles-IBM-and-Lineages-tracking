%%%05-09-24 WORKING VERSION NEW FORWARD TRACKING SCRIPT
%%%PART 1: Identify cell lineages by closest distance and correct for 
%%%cells that find no immediate neighbours 
%%%version for Anthony%%%
% this version starts with the 'final_merged_table.xlsx'
% output is 'tracked_cells_Matlab.csv' and 'tracked_lineages_matlab.png'

% don't forget to go to line 94 to set the time of entry to stationary
% phase
% given the poor quality of stationary phase segmentation, for now this is
% limited to the time points of exponential growth, e.g., 1-t_stat


% Script has four parts 
% Part one: define the most plausible draft next kin
% for each cellID position based on shortest distance Then try to correct
% this by only keeping those next kin which have a follow up in the t+2
% round 
% Part two: build a draft lineage path for each founder cell 
% Save an output file with the lineage and daugther cell ID information and plot a
% summary PNG. 
% Part three: Build a draft generation table and a mother ID
% table for each lineage and try to reiteratively correct this by finding
% branches on the tree that may appear too soon (based on 10th percentile)
% and may be placed by distance on branches that seem too slow (based on
% 90th percentile).
% Part four: make the final generation table, the mother ID table, the
% lineage plot table and tables that can link to GFP and CHE values

% version for a single position
% here is an example

% cd('/Volumes/RECHERCHE/FAC/FBM/DMF/jvanderm/ice-clc/D2c/RawData/Anthony_C/Microscopy/Timelapse/231206_07fpr/_1/Pos1')

%% Section one: look for errors in positions - exclude those positions with time points with a single cell

clear;
close all;

% general parameters

threshold = 25; %distance to chose above which we consider it unlikely to have a neighbours or daughters

% verify if there is a tracking table and read the data

all_pos = readtable(strcat('Data/','final_merged_table_pos14.xlsx'));

%sort all rows by time point and make an array in Matlab

all_pos=sortrows(all_pos,'Timepoint');

% keep only the first four columns with the order 'MotherID','geomX','GeomY','Timepoint'

data = all_pos(:,[1:4, 8:9, 11]);

% transform to array

times = unique(data(:,4));
times = table2array(times);
data = table2array(data);

for j=1:size(times,1)-1

    s(j)=length(find(data(:,4)==times(j))); 

end

%%remove lines with only one cell (which would stop the script);

if ~isempty(find(s==1, 1))

    single_lines = find(s==1);

	for i = 1:length(single_lines)
	    data(data(:,4) == times(single_lines(i)),:) = [];
    end

end

% prepare an updated 'times' variable

times = unique(data(:,4));

% make a growth plot to define the approximate timing of the stationary
% phase which you can report in line 95

figure(1)
plot(smooth(s))

%% collect the cells for tracking and define their most plausible 'next_kin'
% based on minimal distance across two images - and then correction for
% those cells where there are three or more closest distances - further
% correction by comparing cells to the timepoint t+1, in order to see
% whether cells missed in t-1 have a related ID at t+1
% define the second closest distance depending on the strain, between 5 and
% 10% missing cells compared to the founder cell positions and numbers
close all

t_stat = 35;
%define the distance matrices (D) and the closest neighbour matrices (I)
 
D = {}; I = {};

%make new empty table with length of data to fill the IDs of next of kin
%collect the data at each time point to make a new table that includes new cellIDs

next_kin = []; 
new_data = [];

%go through all time points and collect successively data from current
%time=1 (the s-array), time+1 (the t-array) and time +2 (t2) for comparison

for j=1:t_stat %size(times,1)-2

ind_t_1 = find(all_pos.Timepoint == times(j));
ind_t_2 = find(all_pos.Timepoint == times(j + 1));
ind_t_3 = find(all_pos.Timepoint == times(j + 2));
s = data(ind_t_1, :);
t = data(ind_t_2, :);
t2 = data(ind_t_3, :);
all_cellIDs = sort(unique(data(:, 1)));

Pos_X_1 = [s(:,2), s(:,3)]';
Pos_X_2 = [t(:,2), t(:,3)]';
Pos_X_3 = [t2(:,2), t2(:,3)]';
Dist = distEuclid(Pos_X_1, Pos_X_2);
Dist_23 = distEuclid(Pos_X_2, Pos_X_3);
Dist_13 = distEuclid(Pos_X_1, Pos_X_3);

vect_Cell_length_1 = s(:, 6); 
height_cell_1 = s(:, 5);
cell_angle_1 = s(:, 7);
vect_Cell_length_2 = t(:, 6); 
height_cell_2 = t(:, 5);
cell_angle_2 = t(:, 7);
vect_Cell_length_3 = t2(:, 6); 
height_cell_3 = t2(:, 5);
cell_angle_3 = t2(:, 7);

[all_pos, vect_Cell_length_1, height_cell_1, cell_angle_1] = Fill_NaN(ind_t_1, Dist', all_pos, vect_Cell_length_1, height_cell_1, cell_angle_1, Pos_X_1, vect_Cell_length_2, height_cell_2, cell_angle_2);
n_1 = height(s(:, 1));
Seg_tot_1 = arrayfun(@(x) Rect2Seg([Pos_X_1(1,x) Pos_X_1(2,x) vect_Cell_length_1(x) height_cell_1(x)], cell_angle_1(x)),1:n_1,'UniformOutput',false); %Find an alternative
cell_ID_1 = s(:, 1);

[all_pos, vect_Cell_length_2, height_cell_2, cell_angle_2] = Fill_NaN(ind_t_2, Dist, all_pos, vect_Cell_length_2, height_cell_2, cell_angle_2, Pos_X_2, vect_Cell_length_1, height_cell_1, cell_angle_1);
n_2 = height(t(:, 1)); 
Seg_tot_2 = arrayfun(@(x) Rect2Seg([Pos_X_2(1,x) Pos_X_2(2,x) vect_Cell_length_2(x) height_cell_2(x)], cell_angle_2(x)),1:n_2,'UniformOutput',false); %Find an alternative
cell_ID_2 = t(:, 1);

[all_pos, vect_Cell_length_3, height_cell_3, cell_angle_3] = Fill_NaN(ind_t_3, Dist_23, all_pos, vect_Cell_length_3, height_cell_3, cell_angle_3, Pos_X_3, vect_Cell_length_2, height_cell_2, cell_angle_2);
n_3 = height(t2(:, 1));
Seg_tot_3 = arrayfun(@(x) Rect2Seg([Pos_X_3(1,x) Pos_X_3(2,x) vect_Cell_length_3(x) height_cell_3(x)], cell_angle_3(x)),1:n_3,'UniformOutput',false); %Find an alternative
cell_ID_3 = t2(:, 1);

h_mat = H_mat_fun(Dist, n_1, n_2, vect_Cell_length_1, vect_Cell_length_2, 1.2*height_cell_1, 1.2*height_cell_2, Seg_tot_1, Seg_tot_2);
h_mat(h_mat < 0) = 0;

h_mat_23 = H_mat_fun(Dist_13, n_1, n_3, vect_Cell_length_1, vect_Cell_length_3, 1.2*height_cell_1, 1.2*height_cell_3, Seg_tot_1, Seg_tot_3);
h_mat_23(h_mat_23 < 0) = 0;

% D has the individual distances between cells in s and in t, and s and t2

D{j} = distEuclid([s(:,2),s(:,3)]',[t(:,2),t(:,3)]');
D{j+1} = distEuclid([s(:,2),s(:,3)]',[t2(:,2),t2(:,3)]');

% columns in I are the ranks of the cells in t, values are the ranks of the
% cells in s; note that we do not restrict the distance to a specific value
% here but take the smallest one. 

[E{j}, I{j}] = min(D{j}); %max(h_mat);%

%Create a function to see if two daugther have exactly the same mother
%distance.
% We impose the distance restriction, but use this later.
E{j}(E{j} > threshold) = nan; %E{j}(E{j} == 0) = nan; %
[E{j+1}, I{j+1}] =  min(D{j + 1}); %max(h_mat_23); %
E{j+1}(E{j+1} > threshold) = nan; %E{j + 1}(E{j + 1} == 0) = nan; %

% define the arrays for the comparison, and correct a 'cells' array to
% identify cases of multiple closest distances, which can then be corrected
% below to build the next_kin list. The 'tells' array is a help to identify
% the columns of cells in the t-list with multiple closest partners in s.
% The 't2ells' array is a help to identify the columns of cells in the
% t2-list with multiple closest partners in s.

cells = s(I{j});
c2ells = s(I{j+1});
tells = t(:,1);
t2ells = t2(:,1);

%Initialiye with nan
hlp=[]; hlp2=[]; %help matrix to collect the values for the next_kin table

% count the occurrences of the cell IDs in the smallest distance comparison

I_u = unique(I{j}); %unique cell IDs in the smallest distance comparison
I_c = []; %collect the occurrences of each of the cell IDs

for i=1:length(I_u)

    I_c(i)=sum(I{j} == I_u(i));

end

% there may be multiple occasions where there are three or more closest
% distances this is solved in the loop below where we can correct the
% 'cells'-matrix by finding new positions of 'next-closest' neighbouring
% cells that occur at permissive sites in C_c

% count the corresponding occurrences in the t=1 's'-array, to define
% permissive sites for permutation

if ~isempty(find(I_c>2, 1))

	I_ID=[];
	
	%find the positions in the I_c matrix which occur more than twice
	I_ID=I_u(I_c>2); 

	%go through those positions individually
	for zz=1:length(I_ID)


	C_c=[]; %collect the occurrences of each of the cell IDs in the s-array in column 1

	for i=1:length(s)

	    C_c(i)=sum(s(i)==cells);

	end

	% identify the nonpermissive sites in C_c where the counts are either 2 or larger, 
	% because these cannot be used for finding new partners

	C_c=C_c>1; %this will be a logical array, will be 'false' where values are 1 or lower


	%make an index of the s-vector to identify the cell IDs in s
	index=unique(s(:,1),'stable');

	%get the positions in t-array, which occur multiple times
	temp=tells(cells==index(I_ID(zz))); 

	%gets the cellIDs of the cells in t
	wer=find(ismember(tells,temp)); 

	%get the corresponding distance matrix portion - the distance columns
	werd=D{j}(:,wer); 

	%subtract the line with the current lowest distances, which we want to replace
	tmp2=werd-werd(I_ID(zz),:); 

	% this line is now 'zero' - find the line and suppress it
	line_zero=find(all(tmp2 ==0,2));
	tmp2(line_zero,:)=nan;

	% now compare to the permissive matrix of s (C_c) to identify potential
	% new partners we suppress those positions by changing to nan and they
	% cannot be found as minimum later on
	non_permissive=find(C_c==1);
	tmp2(non_permissive,:)=nan;
	
	% find the next two lowest distances, keep row numbers
	% make a copy of tmp2 for comparison
	tmp3=tmp2;
	[tmp2,r_tmp]=mink(reshape(tmp2,[],size(temp,1)),2); 

	% verify if the lowest distances have duplicate rows, 
	% but only when there are 4 or more columns in tmp2

	if size(tmp2,2)>3 && sum(double(ismember(histcounts(r_tmp(1,:),'BinEdge',[unique(r_tmp(1,:)),max(r_tmp(1,:))+1]),2)))>0
	
		%find double ones and collect their row numbers

		[val,edg]=histcounts(r_tmp(1,:),'BinEdge',[unique(r_tmp(1,:)),max(r_tmp(1,:))+1]);

		repeatedElements = edg(val >= 2);

		%in very rare cases, there may be two or more repeatedElements. this
		%cannot be resolved to avoid an error, we only take the first of these
		%(arbitrarily). Isaline works on a different solution here.

		[r,c]= find(r_tmp(1,:) == repeatedElements(1));

		%c are the column numbers of where the repeats occur
		%now interchange those
		r_tmp(r(1),c(2))=r_tmp(2,c(2));
		
		%don't change the value in tmp2, but now delete the second line
		%to make the array suitable for continuation below

		r_tmp(2,:)=[];
		tmp2(2,:)=[];

	else
		r_tmp(2,:)=[];
		tmp2(2,:)=[];
	end
	
	%sort to find the potential new partner cells, but keep sorting order
	[tmp2, order_tmp]=sort(tmp2); 

	% go through each one of those to replace the s-cellID in the 'cells'
	% list subtract by 2 (paired cells that are kept) to know how many values
	% have to be changed
	change = size(temp,1)-2; 
	for kk=1:change
		%coordinates defined by column and row from the sorted order			
		coordinate=[order_tmp(kk),r_tmp(order_tmp(kk))];
		%find the appropriate column
		wer_column=wer(coordinate(:,1));
		%find the new partner cell ID in the s-list
		s_cell=s(coordinate(:,2));
		%replace the cellID in the cells-list
		cells(wer_column)=s_cell;
	end
	end
end


% The corrected 'cells'-matrix is now used to build the next_kin table
% where we give only three possibilities no neighbour detected (both values
% become nan), one neighbour detected (no cell division; [cellID,nan]) or
% two neighbours (means: cells divide in two daughters).
% Values then become the new cellIDs. We get back to the s-index list with
% order maintained. The hlp-list from above is now used to build the
% next_kin list

index=unique(s(:,1),'stable'); %important not to change the order

% now limit the list to those with distances below threshold

tells(isnan(E{j}))=nan;
%Create the hlp matrix
for k=1:length(index)
    if sum(cells==index(k))==1
        hlp(k,:)=[tells(cells==index(k)),nan];	
    elseif sum(cells==index(k))==2
        tmp_row=tells(cells==index(k));
        hlp(k,:)=sort(tmp_row,'ascend');
    else
        hlp(k,:)=[nan,nan];
    end
end


index2=unique(s(:,1),'stable'); %important not to change the order
%only limit the list to those with distances below threshold
t2ells(isnan(E{j+1}))=nan; 
for k=1:length(index2)
    if sum(c2ells==index2(k))==1
        hlp2(k,:)=[t2ells(c2ells==index2(k)),nan];
    elseif sum(c2ells==index2(k))==2
        tmp_row=t2ells(c2ells==index2(k));
        hlp2(k,:)=sort(tmp_row,'ascend');
    else
        hlp2(k,:)=[nan;nan];
    end
end


% pruning the hlp array: can we link s to t and identify cells with double nan that
% are missing at time=t

two_nan_in_t=index(~any(~isnan(hlp), 2));

%can we find the cell ID identities of those that may be within 25 px  - look in D?

pot_daughter_dist_for_2nan_in_t=D{1,j}(~any(~isnan(hlp), 2),:)<threshold;

% and what are the identities of the dividing cells in hlp?

divided_in_t=hlp(~any(isnan(hlp), 2),:);

% what are their current closest distances?

% to get those cellIDs; can be multiple lines

for zz=1:size(pot_daughter_dist_for_2nan_in_t,1)

    % find their cell IDs and if they are part of the dividing cells in t2

    tmp=tells(pot_daughter_dist_for_2nan_in_t(zz,:));

    % are they part of this?

    tmp=tmp(ismember(tmp,divided_in_t));

	    if ~isempty(tmp)

	        % find their distances to the 2nan cell in t2
    
	        tmp_dist=D{1,j}(index==two_nan_in_t(zz),ismember(tells,tmp));
    
	        % what is their current minimal distance
    
	        curr_dist=min(D{1,j}(:,ismember(tells,tmp)));
    
	        % take the difference
    
	        tmp_diff=abs(curr_dist-tmp_dist);
    
	        % find the position of the minimal difference
    
	        tmp_diff=tmp_diff==(min(tmp_diff));
    
	        % Find the corresponding cell to exchange.
            % If by chance both distances are the same, and tmp_diff would
            % have two values, we can arbitrarily take only one, to avoid
            % that the program gets stuck.
    
	        cellID_exchange=tmp(tmp_diff);
            cellID_exchange=cellID_exchange(1);
    
	        % replace the cellID at this position in hlp by a nan and reorder to
	        % have Daughter 2 becoming Daughter1
    
	        [r,c]=find(hlp==cellID_exchange);
    
	        tmp_row=hlp(r,:);
	        tmp_row(c)=nan;
	        hlp(r,:)=sort(tmp_row,'ascend');
    
	        % and add the cellID to the position to change in the previous 2xnan row
    
	        hlp(ismember(index,two_nan_in_t(zz)),1)=cellID_exchange;

	    end

clearvars tmp_row r c cellID_exchange tmp_diff curr_dist tmp_dist

end

clearvars two_nan_in_t pot_daughter_dist_for_2nan_in_t divided_in_t

% Now link s to t+1 and identify cells with double nan that are missing
% at time=t, but may link at t2

two_nan_in_t=index(~any(~isnan(hlp), 2));

%verify if the cell can be linked at t2

two_nan_in_t2=index(~any(~isnan(hlp2), 2));

two_nan_in_t=setdiff(two_nan_in_t,two_nan_in_t2);

if ~isempty(two_nan_in_t)

% make new cell IDs

    newIDs = max(all_cellIDs)+[1:size(two_nan_in_t)]';

    %Remove 
    all_cellIDs = max(newIDs);

    % find those lines in the data file to add them back but with a new cell ID
    
    updata=data(ismember(data(:,1),two_nan_in_t),:);
    updata(:,1)=newIDs;
    
    % don't forget to increase the time
    
    updata(:,4)=updata(:,4)+1;
    
    data=[data;updata];
    
    % then replace the IDs in the hlp file corresponding rows
    
    %find the positions and replace
    
    hlp(ismember(index,two_nan_in_t))=newIDs;

end

%%this defines and builds the next_kin table

next_kin=vertcat(next_kin,hlp);
new_data=vertcat(new_data,s);

clearvars temp t_comp two_nan_in_t newIDs updata 

end

%%%%%%%end of the first section of the cell tracking optimization%%%%%%%%
%data has the order 'MotherID','geomX','GeomY','Timepoint' next_kin table
%has next_ID,nan or daughterID1 and daughterID2
%% Section: Find founder cells.
% Given the difficulty to establish the true founder cells at the first
% image, we look for unique mothers that can make a lineage within the
% first 5 time points

pp=1;
missing_lineage={};

%combined new_data and next_kin

all_tracks=[new_data(:,1:4),next_kin];

%first 5 timepoints

sub=all_tracks(ismember(all_tracks(:,4),times(1:5)),:);

missing_table=sub;

while length(missing_table)>1%0.1*length(sub)

% to collect the lineages
% basically we try to create a path from the mother to all the descendants,
% including those that have no further offspring

missing_founder_cells=missing_table(1,1);
time= missing_table(1,4);

fcID=missing_founder_cells;

%find the corresponding row and column numbers in the data file

[r,c]= find(missing_table(:,1)==fcID);

if ~isempty(r)

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

	motherID=find(sub(:,1)==dcID); 

	%find the new corresponding row numbers

	[r,c]=find(missing_table(:,1)==dcID); 

	%go to the next time point

	time = missing_table(r,4); 

	%add the new IDs to the lineage path

	cellID=[cellID, dcID]; 

	%empty the list of daughter IDs for the next generation

	dcID=[];
end

%collect the lineage for this founder cell and start again
end

missing_lineage{pp}=cellID;
missing_cells=setdiff(missing_table(:,1),[missing_lineage{:}]);
missing_table=missing_table(ismember(missing_table(:,1),missing_cells),:);

pp=pp+1; %add to the tracker

end

%take from every new lineage the first motherID, which forms the founder cells

founder_cells=[]; 

for kk=1:length(missing_lineage)

founder_cells(kk)=missing_lineage{kk}(1);

end


%%  Section 3: Find draft lineages %%%%%%%%%%%

% we combine the updated data and the next kin in a new array

all_tracks=[new_data(:,1:4), next_kin];

% find lineages
% loop through all cells, lineage by lineage until the table has become
% smaller than 10% of its original length. Hence: 90% of cells are tracked

pp=1;
missing_lineage={};

missing_table=all_tracks;

while length(missing_table) > 0.1*length(next_kin)

missing_founder_cells = missing_table(1,1);
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

	cellID=[cellID, dcID]; 

	%empty the list of daughter IDs for the next generation

	dcID=[];
end

%collect the lineage for this founder cell and start again

missing_lineage{pp}=cellID;
missing_cells=setdiff(missing_table(1:length(missing_table),1),[missing_lineage{:}]);
missing_table=missing_table(ismember(missing_table(:,1),missing_cells),:);

pp=pp+1; %add to the tracker

end

%%list lineages in an array, every column is a lineage

lng=[]; %to collect the lineages in the array

for kk = 1:length(missing_lineage)

%give every new lineage a number that corresponds to the next subfile in
%the lineage cell array by using the kk-index

lng(:,kk) = double(ismember(all_tracks(:,1),missing_lineage{kk}))*kk;

end

% combine all into a single lineage column where every number is a unique
% lineage

lng_tot = sum(lng,2);


% add to the all_tracks file as the 7th column

all_tracks(:,7)=lng_tot(1:length(next_kin));

% this is the output of this part and saved as .csv file


Variables={'MotherID','GeomX','GeomY','Timepoint','DaughterID1','DaughterID2','Lineage'};
t=array2table(all_tracks,'VariableNames',Variables);
writetable(t,'tracked_cells_Matlab.csv');


%% Section: make an intermediate plot of the tracked lineages
close all

figH=figure;

subplot('Position',[0.1,0.1,0.5,0.5]);

lng_nrs=unique(all_tracks(:,7)); %all lineage numbers

for i=1:length(lng_nrs)

tmp=all_tracks(all_tracks(:,7)==lng_nrs(i),:);

scatter(tmp(:,2),tmp(:,3),'.')
hold on

end

saveas(figH,'tracked_lineages_matlab.png','png');
