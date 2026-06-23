%% Part four: connect stationary phase data and make final tables
% In the first section we try to connect stationary phase data 
% we use a simplified version of the distance matrix, keeping the 
% last timepoint of the previous tracking as a mask for every other time point

% In the second section we combine the stationary phase tracks
% to the rebuild mother ID table.
% This table is then used to find the corresponding GFP and mCHE
% values based on the corresponding cellID.
% We also make a plot for the genealogy tree and calculate
% the generations for each of the cell IDs in the same format.

% We start from the Lineages and the cellcounts structures
% You can open them from the saved version or continue from part 3.

clearvars -except cellcounts Lineages

% We also open the Dimalis file with the fluorescence values to compare
% with. 

features = readtable('final_merged_table.xlsx');

% data has the order 'MotherID','geomX','GeomY','Timepoint'

data=sortrows(features,'Timepoint');

times=unique(data(:,4));
times=table2array(times);
data=table2array(data);

% go through the individual lineages

for i=1:length(Lineages)

% define the stat phase mask of the cells

stat_mask=Lineages{1,i}(Lineages{1,i}.Timepoint==max(Lineages{1,i}.Timepoint),:);
stat_tp=max(Lineages{1,i}.Timepoint);

%define the distance matrices (D) and the closest neighbour matrices (I)
 
D={};I={};

%make new empty table with length of data to fill the IDs of next of kin
%collect the data at each time point to make a new table that includes new cellIDs

stat_ID_table=stat_mask.MotherID;

stat_mask_geo=[stat_mask.GeomX,stat_mask.GeomY];


%%go through all stationary phase time points and collect successively time=1, time+1 for distance comparison

for j=stat_tp:size(times,1)-2

t=[]; %initiate fresh arrays

t=data(data(:,4)==times(j+1),:);

%D has the individual distances between cells in s and in t

D{j} = pdist2([t(:,2),t(:,3)],stat_mask_geo);
[~, I{j}] = pdist2([t(:,2),t(:,3)],stat_mask_geo,'euclidean','Smallest',1);

%columns in I are the IDs of the cells in stat_mask
%the ones that are too far away become nan

I{j}(min(D{1,j},[],1)>25)=nan;

% identify the new cell IDs- can't use logic indexing because sometimes nan

cellIDs=[];

for zz=1:length(I{j})

	if ~isnan(I{j}(zz))

	cellIDs(zz)=t(I{j}(zz));

	else
	
	cellIDs(zz)=nan;

	end

end

stat_ID_table=[stat_ID_table,cellIDs'];

end

%save into the cellcounts structure

cellcounts.lineage(i).stat_ID_table=stat_ID_table;

end


clearvars -except Lineages cellcounts features stat_tp
%% Remove lineages that are too small
% and that block the script later on. Here the limit is set to 50 IDs (line
% 105)
load('Lineages.mat'); load('cellcounts.mat');

for i=1:size(Lineages,2)

	if size(Lineages{i},1)<50
	    Lineages{i}=[];
	    cellcounts.lineage(i).gen_table=[];
	    cellcounts.lineage(i).stat_ID_table=[];
	end
end
%% Section two
% go through all Lineages to make the motherID table, a generation
% time table, a plot table, and the corresponding GFP and mCHE tables

stat_tp=max(Lineages{1,i}.Timepoint);
features = readtable('final_merged_table.xlsx');
for z=1:size(Lineages,2)

sub=Lineages{z};

if ~isempty(sub)

% correct MotherID for Mother_ID

sub=renamevars(sub,'MotherID','Mother_ID');

% make a Mother ID table
% convert generation table to time by all cell IDs,
% shape by time_length

%time_length=(1:max(sub.Timepoint));
time_length=(1:stat_tp);

mother_ID_table = zeros(size(cellcounts.lineage(z).gen_table,1),size(time_length,2));
mother_ID_table(mother_ID_table == 0) = nan;

% now fill this successively with each mother ID and its path to the next division
% for this we need the sub table with all intermediate cell IDs

lin=cellcounts.lineage(z).gen_table;

for i=1:size(lin,1)

	ID_birth=[]; ID_death=[]; dID=[]; r=[]; r2=[]; aa=[]; hlp=[];  tmp=[];
    cell_birth_tp=[]; cell_death_tp=[];
	ID_birth=lin.ID_birth(i);
	ID_death=lin.ID_death(i);

    %This loop is totally useless. It is exactly the same condition twice.
	if lin.cell_birth_tp==0 %make this artifically time 1

	mother_ID_table(i,lin.cell_birth_tp(i)+1)=ID_birth;

	else
	
	mother_ID_table(i,lin.cell_birth_tp(i)+1)=ID_birth;

	end

	%now see if ID_birth and ID_death are the same, which means that the first cell
	%has divided at this point and has disappeared. Then we only want the 

	if ID_birth==ID_death

	dID=ID_birth;

	else

	dID=sub.DaughterID1(sub.Mother_ID==ID_birth);

	end

	while dID~=ID_death %

	[r,c]=find(sub.Mother_ID==dID);
	mother_ID_table(i,sub.Timepoint(sub.DaughterID1==dID)+1)=dID;

		%if happens that the final daughter ID is the same as what we are looking for
		if dID==sub.DaughterID1(r)
			break
		else	
		dID = sub.DaughterID1(r);
		end
	end

	if ~isnan(ID_death)

	mother_ID_table(i,lin.cell_death_tp(i)+1)=ID_death;

	else %stop here because the cell has disappeared

	mother_ID_table(i,sub.Timepoint(sub.Mother_ID==ID_birth)+1)=ID_birth;
	end

end


% find the mother ID that belongs to the ID_birth in the current search
% and work the table backwards to fill the rows with the previous generations
% go through all the rows of the mother_ID_table, except the first, which is complete already

% remove rows with all nan

mother_ID_table(~any(~isnan(mother_ID_table), 2),:)=[];

%remove values that are zero

mother_ID_table(mother_ID_table==0)=nan;

%fill all the lineages with cellIDs backwards

for i=2:size(mother_ID_table,1)

	tmp=mother_ID_table(i,~isnan(mother_ID_table(i,:)));
	start_ID=tmp(1); %the first ID before the new cell lineage
	[r,c]=find(mother_ID_table(i,:)==start_ID);

	lineage_ID=sub.Mother_ID(find(sub.DaughterID1 == start_ID | sub.DaughterID2== start_ID)); 

	if ~isempty(lineage_ID>1)

		lineage_ID=lineage_ID(1);

	end

	mother_ID_table(i,c(1)-1)=lineage_ID;

	%identify again the first cell ID in the table as a reference to fill

	tmp_start=min(mother_ID_table(1,~isnan(mother_ID_table(1,:))));

	while lineage_ID~=tmp_start

	%place this cell ID at the appropriate position before the starting hlp_ID

	[r2,c2]=find(mother_ID_table(i,:)==lineage_ID);
	
	lineage_ID=sub.Mother_ID(find(sub.DaughterID1 == lineage_ID | sub.DaughterID2== lineage_ID)); 
	mother_ID_table(i,c2-1)=lineage_ID;

	end
end

% now we have to fuse the stat phase ID table to the proper positions in the 
% Mother ID table

% Extend the mother_ID_table with the width of the stat table minus 2
% because the first two columns in stat table are a repetition of the last of
% the mother_ID time points

hlp=zeros(size(mother_ID_table,1),size(cellcounts.lineage(z).stat_ID_table,2)-2);
hlp(hlp==0)=nan;
mother_ID_table=[mother_ID_table,hlp];

% find the rows in the mother_ID_table with the last cell ID being the same
% as the first in the stat array

rows=[]; cols=[];

[rows,cols]=find(ismember(mother_ID_table(:,stat_tp:41),cellcounts.lineage(z).stat_ID_table(:,1)));
rows=unique(rows); %to prevent accidental doubles

% identify corresponding rows and complete each of those rows in the
% mother_ID_table with the cell IDs from the stat phase

for ii=1:length(rows)

mother_ID_table(rows(ii),(stat_tp+cols(ii)):(end-1))=cellcounts.lineage(z).stat_ID_table(cellcounts.lineage(z).stat_ID_table(:,1)==mother_ID_table(rows(ii),stat_tp+cols(ii)-1),3:2+length(mother_ID_table(rows(ii),(stat_tp+cols(ii)):(end-1))));

end 

	
% add to the cellcounts structure

cellcounts.lineage(z).mother_ID_table=mother_ID_table;

% now find the unique mother IDs in the table and combine to the relevant features
% to fill the GFP and CHE tables

tmp=features(ismember(features.Mask_nb,unique(mother_ID_table(~isnan(mother_ID_table)))),:);

% make similar eGFP table

GFP_table=zeros(size(mother_ID_table,1),size(mother_ID_table,2));
GFP_table(GFP_table==0)=nan;

for i=1:size(mother_ID_table,1)

	for j=1:size(mother_ID_table,2)

	X=mother_ID_table(i,j);
	if ~isnan(X)
	hlp=tmp.GFP_mean(tmp.Mask_nb==X);
		if ~isempty(hlp)
		GFP_table(i,j)=hlp;
		else
		GFP_table(i,j)=GFP_table(i,j-1);
		end
	else
	GFP_table(i,j)=nan;
	end

	end
end

% add to the cellcounts structure

cellcounts.lineage(z).GFP_table=GFP_table;

% make similar mCHE table

CHE_table=zeros(size(mother_ID_table,1),size(mother_ID_table,2));
CHE_table(CHE_table==0)=nan;

for i=1:size(mother_ID_table,1)

	for j=1:size(mother_ID_table,2)

	X=mother_ID_table(i,j);
	if ~isnan(X)
	hlp=tmp.mCherry_mean(tmp.Mask_nb==mother_ID_table(i,j));
		if ~isempty(hlp)
		CHE_table(i,j)=hlp;
		else
		CHE_table(i,j)=CHE_table(i,j-1);
		end
	else
	CHE_table(i,j)=nan;
	end

	end
end

% add to the cellcounts structure

cellcounts.lineage(z).CHE_table=CHE_table;


% make similar generation table

gen_table=zeros(size(mother_ID_table,1),size(mother_ID_table,2));
gen_table(gen_table==0)=nan;

for i=1:size(mother_ID_table,1)

	for j=1:size(mother_ID_table,2)

	X=mother_ID_table(i,j);
	if ~isnan(X)
	hlp=lin.generation(find(lin.ID_birth==mother_ID_table(i,j)));
		if ~isempty(hlp)
		gen_table(i,j)=hlp;
		elseif X==mother_ID_table(1,1) %exception for first cell
		gen_table(i,j)=1;
		else
		gen_table(i,j)=gen_table(i,j-1);
		end
	else
	gen_table(i,j)=nan;
	end

	end
end

% add to the cellcounts structure

cellcounts.lineage(z).generation_table=gen_table;

% make similar generation plot table
% split into two parts, first until stat phase
% then from stat phase to end

plot_table=zeros(size(mother_ID_table,1),size(mother_ID_table,2));
plot_table(plot_table==0)=nan;

for i=1:size(mother_ID_table,1)
    if size(plot_table,1)==1
        break
    end
	for j=1:stat_tp%

	X=mother_ID_table(i,j);
	if ~isnan(X)
	hlp=lin.generation(find(lin.ID_birth==mother_ID_table(i,j)));
		if X==tmp_start %first mother cell
		plot_table(i,j)=2^0; %y-pos starting point = 1
	
		elseif ~isempty(hlp) && X~=tmp_start && X~=mother_ID_table(i,j-1)
			
			if ismember(X,lin.ID_birth) && ismember(X, sub.DaughterID1)
				plot_table(i,j)=plot_table(i,j-1)+2^(0-hlp);
			elseif ismember(X, lin.ID_birth) && ismember(X,sub.DaughterID2)
				plot_table(i,j)=plot_table(i,j-1)-2^(0-hlp);
			else
				plot_table(i,j)=plot_table(i,find(~isnan(plot_table(i,:)),1,'last'));
			end

        	else
			plot_table(i,j)=plot_table(i,j-1);
		end
	else
	plot_table(i,j)=nan;
	end

    end
end

for i=1:size(mother_ID_table,1)
    if size(plot_table,1)==1
        break
    end
	for j=stat_tp:size(mother_ID_table,2)

	X=mother_ID_table(i,j);
	    if ~isnan(X)
	    hlp=lin.generation(find(lin.ID_birth==mother_ID_table(i,j)));
		    if X==tmp_start %first mother cell
		    plot_table(i,j)=2^0; %y-pos starting point = 1
	
		    elseif ~isempty(hlp) && X~=tmp_start && X~=mother_ID_table(i,j-1)
			
			    if ismember(X,lin.ID_birth) && ismember(X, sub.DaughterID1)
				plot_table(i,j)=plot_table(i,j-1)+2^(0-hlp);
			    elseif ismember(X, lin.ID_birth) && ismember(X,sub.DaughterID2)
				plot_table(i,j)=plot_table(i,j-1)-2^(0-hlp);
			    else
				plot_table(i,j)=plot_table(i,find(~isnan(plot_table(i,:)),1,'last'));
			    end

        	else
			plot_table(i,j)=plot_table(i,j-1);
		    end
	    else
	    plot_table(i,j)=plot_table(i,j-1);
	    end

    end
end

% add to the cellcounts structure

cellcounts.lineage(z).plot_table=plot_table;
end
end

save Lineages.mat Lineages
save cellcounts.mat cellcounts

%% example plot of the lineage tree

z=1; %chose the lineage to plot

close all

figure(1)
plot(1:size(cellcounts.lineage(z).plot_table,2),cellcounts.lineage(z).plot_table,'k');

figure(2)

plot(1:size(cellcounts.lineage(z).CHE_table,2),cellcounts.lineage(z).CHE_table)

figure(3) %qqplot at time 50

qqplot(cellcounts.lineage(z).CHE_table(:,50))
