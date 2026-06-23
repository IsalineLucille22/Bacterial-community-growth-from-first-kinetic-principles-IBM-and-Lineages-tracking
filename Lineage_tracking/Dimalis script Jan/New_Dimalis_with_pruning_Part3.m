%%%05-09-24 WORKING VERSION NEW FORWARD TRACKING SCRIPT version for
% Anthony
%%%%% PART 3: make a draft gen_table and mother_ID table to start the
% pruning
% This part starts with the file 'stitched_tracked_cells_Matlab.csv'. It
% will maintain two Matlab outputs: the Lineages and the cellcounts
% structure. The cellcounts will get the final tables in part 4 for
% plotting.

% Important: define the frequency of the images taken in line 41

%% Section one: read the lineages and make a draft gen_table and mother_ID_table
 
% collect all the lineages

clear
close all

stitching=readtable('stitched_tracked_cells_Matlab.csv');

lng=unique(stitching.Lineage);

Lineages={};

for i=1:length(lng)
Lineages{i}=stitching(stitching.Lineage==lng(i),:);
end

clearvars -except Lineages

% define distance threshold to reconnect branches

thr = 30; %in pixel distance

%% Section: go through all the lineages individually and make draft gen_tables
% The draft gen_table is then updated every round
% depending on the branches of the tree that can be replaced elsewhere.

for k=1:size(Lineages,2)

sub=Lineages{k};

% what is the time frequency (in per hour) of the images

time_freq=3;

st=sub.MotherID(1); %st is the founder cell from the lineage

gen_table=[]; %to collect the data for the generation table
gen=1; %generation counter

% Go through the whole table and stop when st becomes empty. This is the
% sign that all cells have been attributed

while ~isempty(st)

% to collect the daughter data from the mother cell(s) for the next round of searches

coll=[]; 

% start the loop to find new daughters from the defined mothers in this round

    for j=1:length(st)

    da=[]; %individual daughter(s) from each single mother
    ro=[]; cell_birth=[]; cell_death=[];
    ID_birth=st(j); ID_death=[];

    % r is the row number of the mother cell
    r=find(sub.MotherID==st(j));

    % to collect the connecting rows from mother to its divided daughters, 
    % which we need to define the point of cell death in the end

    cell_rows=[]; 

    % verify if a connecting cell exists and then connect
    % if it doesn't - then break the loop and add only the cell birth values

	    if ~isempty(r)

		%record cell birth as where this cell first appeared in the daughter columns
		[ro,co]=find(sub.DaughterID1==st(j) | sub.DaughterID2==st(j));

			if isempty(ro)
				%this situation arises only for the first founder cell
				cell_birth=sub.Timepoint(r);  
			else
				%sometimes happens that cell is found twice by error
				%then take the earliest timepoint
				if length(ro)>1
					cell_birth=sub.Timepoint(min(ro));
				else
					cell_birth=sub.Timepoint(ro);
				end
			end

		da=[sub.DaughterID1(r),sub.DaughterID2(r)]; 


        %now there are only four possibilities: either the cells divide; in
        %that case there will be two da-values or there is only a single
        %connecting cell: no division yet then the value of the second
        %daughterID is NaN or there are no connecting cells then both
        %values are NaN or there is no new row to be found - end of table
        %reached

		    if ~isnan(da(2)) %immediate new daughter cells
			    cell_rows=[r];

			    %take the same timepoint where this cell divided as its cell_death
			    cell_death=sub.Timepoint(cell_rows); 

			    %take the ID of the same timepoint before the division
			    %ID_death=sub.DaughterID1(cell_rows);  
			    ID_death=sub.MotherID(cell_rows);

			    %collect the values in a temporary array
			    tmp(j,:)=[gen,ID_birth,ID_death,cell_birth,cell_death];


		    elseif ~isnan(da(1)) %now we have to continue to search
			    while isnan(da(2))

			    r=find(sub.MotherID==da(1)); %the row number of the connecting cell
			    da=[sub.DaughterID1(r),sub.DaughterID2(r)]; 

			    %the new value(s) of the daughter(s)
			    cell_rows=[cell_rows;r];

				%if the next round finds no more mothers
				%because end of table, stop

				    if isempty(r)	
					    break
				    end

				%sometimes there can be an error in the list when the mother has the
				%same ID as its daughter this would cause an endless loop, so this has
				%to be detected before continuing

				    if da(1)==sub.MotherID(r)
					    da(1)=NaN;
				    elseif da(2)==sub.MotherID(r)
					    da(2)=NaN;	
				    end


				%if the next round has two NaN, stop

				    if isempty(find(sub.MotherID==da(1), 1))
					    break
				    end


			    end %end of the while loop

			%another check for the case where the immediate next cell will have no
			%daughters and the cell_rows would only be length 1. Then we cannot
			%establish the cell death. We add the number one artifically, that can be
			%subtracted in the following line

			    if isempty(cell_rows)

			    %take the last timepoint of the table as its cell_death
			    cell_death=max(sub.Timepoint); 

			    %take the ID of the last timepoint before the division
			    ID_death=sub.DaughterID1(find(sub.MotherID==st(j)));  

			    %collect the data in the tmp array
			    tmp(j,:)=[gen,ID_birth,ID_death,cell_birth,cell_death];

			%
			    elseif length(cell_rows)==1

			    cell_rows=[cell_rows;1];
	
			    %take the ID of the last timepoint before the division
			    ID_death=sub.MotherID(cell_rows(end-1));  

			    %take the last timepoint before this cell divided as its cell_death
			    cell_death=min(sub.Timepoint(find(sub.DaughterID1==ID_death))); 

			    %collect the data in the tmp array
			    tmp(j,:)=[gen,ID_birth,ID_death,cell_birth,cell_death];

			    %
			    elseif isnan(da(1)) && isnan(da(2))

				    if ~isempty(cell_rows)

				    cell_rows=cell_rows;
				    %take the last timepoint before this cell divided as its cell_death
				    cell_death=sub.Timepoint(cell_rows(end-1)); 

				    %take the ID of the last timepoint before the division
				    ID_death=sub.DaughterID1(cell_rows(end-1));  

				    %collect the data in the tmp array
				    tmp(j,:)=[gen,ID_birth,ID_death,cell_birth,cell_death];


				    else
				    cell_rows=r;

				    %take the same timepoint where this cell divided as its cell_death
				    cell_death=sub.Timepoint(cell_rows); 

				    %take the ID of the same timepoint before the division
				    ID_death=sub.DaughterID1(cell_rows);  

				    %collect the values in a temporary array
				    tmp(j,:)=[gen,ID_birth,NaN,cell_birth,NaN];

				    end

			    %
			    else

			    %take the last timepoint before this cell divided as its cell_death
			    cell_death=sub.Timepoint(cell_rows(end-1)); 

			    %take the ID of the last timepoint before the division
			    ID_death=sub.DaughterID1(cell_rows(end-1));  

			    %collect the data in the tmp array
			    tmp(j,:)=[gen,ID_birth,ID_death,cell_birth,cell_death];

			    end

		else %both values are NaN

			cell_rows=[r];

			%take the same timepoint where this cell divided as its cell_death
			cell_death=sub.Timepoint(cell_rows); 

			%take the ID of the same timepoint before the division
			ID_death=sub.DaughterID1(cell_rows);  

			%collect the values in a temporary array
			tmp(j,:)=[gen,ID_birth,NaN,cell_birth,NaN];

		
		end

	else 
	
	%no connecting cell has been found, we assume it died or disappeared
	%its ID_death and time of cell_death will become NaN
	% cell_birth is where it appeared

	tmp(j,:)=[gen,ID_birth,NaN,sub.Timepoint(sub.DaughterID1==ID_birth | sub.DaughterID2==ID_birth),NaN];
	
	end

	%collecting the newly found daugthers in the loop for the next round

	coll=[coll;da];

%end of the loop for search of the daughters from the previously identified mothers

end

%combine the collected values into the final table
	
gen_table=[gen_table;tmp];

%liberate the temporary array for the next round
tmp=[];

%redefine the newly found offspring into the search for the next round

st=reshape(coll,1,[]);

%remove cells that have no further offspring
st(isnan(st))=[]; 

gen=gen+1; %increase the generation counter

end %end of first while loop to go through all cells until no more mothers

% then add the generation times per cell in a new column (6)

for jj=1:size(gen_table,1)

	if ~isnan(gen_table(jj,5))
	gen_time=(gen_table(jj,5)-gen_table(jj,4))/time_freq; %in hours
	gen_table(jj,6)=gen_time; %in hours

	else
	gen_table(jj,6)=nan;

	end

end


%%collect in the final structure

Variables={'generation', 'ID_birth', 'ID_death', 'cell_birth_tp', 'cell_death_tp','generation_time'};
t=array2table(gen_table,'VariableNames',Variables);

cellcounts.lineage(k).gen_table=t;


clearvars sub gen_table gen
 
end %end of lineages loop (k)

%% Section: now take the draft gen_tables and update the lineages
% Try to improve the mother_ID_table by finding branches that appear too soon
% and compare them to branches that seem too long, to see if they can fit
% by short distance to an existing cell that did not divide and may have been
% misclassified.
% Too soon and too slow is defined from the distribution of the observed
% generation times in the draft table - typically the 10th and the 90th
% percentiles, or (here) for low the value of 0.4. See lines 406-408.
% Finally, it will update the lineages.

for k=1:size(Lineages,2)

    sub=Lineages{k};
    sub=renamevars(sub,'MotherID','Mother_ID');

    % we make a copy of the gen_table to improve

    gen_table=cellcounts.lineage(k).gen_table;

    % make a mother_ID_table that has time as x and generation as y, and is filled
    % by row with connecting cell IDs, starting when they appear

    time_length=(1:max(sub.Timepoint));

    %%what is the time frequency (in per hour) of the images

    time_freq=3;

    mother_ID_table=zeros(size(gen_table,1),size(time_length,2));
    mother_ID_table(mother_ID_table==0)=nan;

    % now fill this successively with each mother ID and its path to the next division
    % for this we need the sub table with all intermediate cell IDs

    lin=cellcounts.lineage(k).gen_table;

    for i=1:size(lin,1)

	    ID_birth=[]; ID_death=[]; dID=[]; r=[]; r2=[]; aa=[]; hlp=[];  tmp=[];
        cell_birth_tp=[]; cell_death_tp=[];
	    ID_birth=lin.ID_birth(i);
	    ID_death=lin.ID_death(i);

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

%remove rows with all nan

mother_ID_table(~any(~isnan(mother_ID_table), 2),:)=[];

%remove values that are zero

mother_ID_table(mother_ID_table==0)=nan;

% this gives the draft mother_ID_table

%look for the generation time distributions to define up and low

nanmedian(gen_table.generation_time);

low_bound=0.4; %for P. putida given its fast division.
%low_bound=prctile(gen_table.generation_time,10); % 10th percentile
hi_bound=prctile(gen_table.generation_time,90); %90th percentile


% what are the cells that divide too fast?

too_fast=gen_table.ID_death(find(gen_table.generation_time<=low_bound));


% Go through each of the too fast cells in the lineage, find the daughters
% and compare whether one of these can be connected elsewhere. Then update
% the gen_table and the mother_ID_table for the next round.

    for j=1:length(too_fast)

	    too_fast_mother=too_fast(j);
	    too_slow=gen_table.ID_birth(find(gen_table.generation_time>=hi_bound));

        too_fast_daughters=[Lineages{k}.DaughterID1(find(Lineages{k}.MotherID==too_fast_mother)),Lineages{k}.DaughterID2(find(Lineages{k}.MotherID==too_fast_mother))];

	    if ~any(~isnan(too_fast_daughters))
		        break

	    else

            tp_too_fast=Lineages{k}.Timepoint(find(Lineages{k}.MotherID==too_fast_mother));

                if tp_too_fast==0 %otherwise gives an error
                    break
                end 

        % use the provisory mother_ID_table to find the cell path of the slow lineage
        % to identify the cell ID at the time of birth of the too_fast daughters

            slow_column=mother_ID_table(:,tp_too_fast); %gives the column
            [slow_lineages,c]=find(ismember(mother_ID_table,too_slow)); %gives the corresponding rows and starting columns

            slow_column=slow_column(slow_lineages); %to combine
            slow_column=slow_column(~isnan(slow_column)); %only show the values that are not nan

	        if ~isempty(slow_column)

            %now calculate the paired distances of the candidates

                slow_column_geo=[Lineages{k}.GeomX(ismember(Lineages{k}.MotherID,slow_column)),Lineages{k}.GeomY(ismember(Lineages{k}.MotherID,slow_column))];
                too_fast_geo=[Lineages{k}.GeomX(ismember(Lineages{k}.MotherID,too_fast_daughters)),Lineages{k}.GeomY(ismember(Lineages{k}.MotherID,too_fast_daughters))];
                if isempty(too_fast_geo)
                    break
                end
                D=pdist2(slow_column_geo,too_fast_geo);

                [r,c]=find(D==min(min(D))); %in the second dimension = rows

                % check if this distance is within the set threshold

	        if min(min(D))>thr %again the distance threshold
		        break
	        else

	            %identify this new mother - make sure there is only one

	            new_mother_ID=slow_column(r(1));

	            %find the corresponding cell to exchange

	            cellID_exchange=too_fast_daughters(c);
                cellID_exchange=cellID_exchange(1,1);

	            % Replace the cellID at its original position by nan. Check
                % that it is not empty or else the script stops

                if ~isempty(cellID_exchange)
	            too_fast_daughters(too_fast_daughters==cellID_exchange)=nan;
	            too_fast_daughters=sort(too_fast_daughters,'ascend');
	            Lineages{k}.DaughterID1(find(Lineages{k}.MotherID==too_fast_mother))=too_fast_daughters(:,1);
	            Lineages{k}.DaughterID2(find(Lineages{k}.MotherID==too_fast_mother))=too_fast_daughters(:,2);

	            % Replace the cellID at the new position in the lineage by a
	            % nan and reorder to have Daughter 2 becoming Daughter1

	            Lineages{k}.DaughterID2(find(Lineages{k}.MotherID==new_mother_ID))=cellID_exchange;
                end
	        end %if loop

	    end %second if loop on the slow column

	end %first if loop

	clearvars too_fast_mother too_slow too_fast_daughters too_fast_geo slow_column slow_column_geo slow_lineages D r c cellID_exchange tp_too_fast

% Reiterate the draft gen_table

sub=[];

sub=Lineages{k};

st=sub.MotherID(1); %st is start

gen_table=[]; %to collect the data for the generation table
gen=1; %generation counter

%go through the whole table and stop when st becomes empty
%this is the sign that all cells have been attributed

while ~isempty(st)

    %to collect the daughter data from the mother cell(s) for the next round of searches

    coll=[]; 

    %start the loop to find new daughters from the defined mothers in this round

    for j=1:length(st)

        da=[]; %individual daughter(s) from each single mother
        ro=[]; cell_birth=[]; cell_death=[];
        ID_birth=st(j); ID_death=[];

        %r is the row number of the mother cell
        r=find(sub.MotherID==st(j));

        %to collect the connecting rows from mother to its divided daughters, 
        %which we need to define the point of cell death in the end

        cell_rows=[]; 

        % verify if a connecting cell exists and then connect
        % if it doesn't - then break the loop and add only the cell birth values

	    if ~isempty(r)

		    %record cell birth as where this cell first appeared in the daughter columns
		    [ro,co]=find(sub.DaughterID1==st(j) | sub.DaughterID2==st(j));

			if isempty(ro)
				%this situation arises only for the first founder cell
				cell_birth=sub.Timepoint(r);  
			else
				%sometimes happens that cell is found twice by error
				%then take the earliest timepoint
				if length(ro)>1
					cell_birth=sub.Timepoint(min(ro));
				else
					cell_birth=sub.Timepoint(ro);
				end
			end

		da=[sub.DaughterID1(r),sub.DaughterID2(r)]; 


        %now there are only four possibilities: 
        %either the cells divide; in that case there will be two da-values
        %or there is only a single connecting cell: no division yet
        %then the value of the second daughterID is NaN
        %or there are no connecting cells
        %then both values are NaN
        %or there is no new row to be found - end of table reached

		    if ~isnan(da(2)) %immediate new daughter cells
			    cell_rows=[r];

			    %take the same timepoint where this cell divided as its cell_death
			    cell_death=sub.Timepoint(cell_rows); 

			    %take the ID of the same timepoint before the division
			    
			    ID_death=sub.MotherID(cell_rows);

			    %collect the values in a temporary array
			    tmp(j,:)=[gen,ID_birth,ID_death,cell_birth,cell_death];


		    elseif ~isnan(da(1)) %now we have to continue to search
			    while isnan(da(2))

			        r=find(sub.MotherID==da(1)); %the row number of the connecting cell
			        da=[sub.DaughterID1(r),sub.DaughterID2(r)]; 

			        %the new value(s) of the daughter(s)
			        cell_rows=[cell_rows;r];

				    %if the next round finds no more mothers
				    %because end of table, stop

				        if isempty(r)	
					        break
				        end

				    %sometimes there can be an error in the list
				    %when the mother has the same ID as its daughter
				    %this would cause an endless loop, 
				    %so this has to be detected before continuing

				        if da(1)==sub.MotherID(r)
					        da(1)=NaN;
				        elseif da(2)==sub.MotherID(r)
					        da(2)=NaN;	
				        end


				    %if the next round has two NaN, stop

				        if isempty(find(sub.MotherID==da(1), 1))
					        break
				        end


			    end %end of the while loop

			%another check for the case where the immediate next cell will
			%have no daughters and the cell_rows would only be length 1
			%then we cannot establish the cell death
			%we add one number artifically, that can be subtracted in the
			%following line

			    if isempty(cell_rows)

			        %take the last timepoint of the table as its cell_death
			        cell_death=max(sub.Timepoint); 

			        %take the ID of the last timepoint before the division
			        ID_death=sub.DaughterID1(find(sub.MotherID==st(j)));  

			        %collect the data in the tmp array
			        tmp(j,:)=[gen,ID_birth,ID_death,cell_birth,cell_death];

			%
			    elseif length(cell_rows)==1

			        cell_rows=[cell_rows;1];
	
			        %take the ID of the last timepoint before the division
			        ID_death=sub.MotherID(cell_rows(end-1));  

			        %take the last timepoint before this cell divided as its cell_death
			        cell_death=min(sub.Timepoint(find(sub.DaughterID1==ID_death))); 

			        %collect the data in the tmp array
			        tmp(j,:)=[gen,ID_birth,ID_death,cell_birth,cell_death];

			%
			    elseif isnan(da(1)) && isnan(da(2))

				    if ~isempty(cell_rows)

				        cell_rows=cell_rows;
				        %take the last timepoint before this cell divided as its cell_death
				        cell_death=sub.Timepoint(cell_rows(end-1)); 

				        %take the ID of the last timepoint before the division
				        ID_death=sub.DaughterID1(cell_rows(end-1));  

				        %collect the data in the tmp array
				        tmp(j,:)=[gen,ID_birth,ID_death,cell_birth,cell_death];


				    else
				        cell_rows=r;

				        %take the same timepoint where this cell divided as its cell_death
				        cell_death=sub.Timepoint(cell_rows); 

				        %take the ID of the same timepoint before the division
				        ID_death=sub.DaughterID1(cell_rows);  

				        %collect the values in a temporary array
				        tmp(j,:)=[gen,ID_birth,NaN,cell_birth,NaN];

                    end

			else

			    %take the last timepoint before this cell divided as its cell_death
			    cell_death=sub.Timepoint(cell_rows(end-1)); 

			    %take the ID of the last timepoint before the division
			    ID_death=sub.DaughterID1(cell_rows(end-1));  

			    %collect the data in the tmp array
			    tmp(j,:)=[gen,ID_birth,ID_death,cell_birth,cell_death];

			end

		else %both values are NaN

			cell_rows=[r];

			%take the same timepoint where this cell divided as its cell_death
			cell_death=sub.Timepoint(cell_rows); 

			%take the ID of the same timepoint before the division
			ID_death=sub.DaughterID1(cell_rows);  

			%collect the values in a temporary array
			tmp(j,:)=[gen,ID_birth,NaN,cell_birth,NaN];

		
		end

	else 
	
	%no connecting cell has been found, we assume it died or disappeared
	%its ID_death and time of cell_death will become NaN
	% cell_birth is where it appeared

	tmp(j,:)=[gen,ID_birth,NaN,sub.Timepoint(sub.DaughterID1==ID_birth | sub.DaughterID2==ID_birth),NaN];
	
	end

	%collecting the newly found daugthers in the loop for the next round

	coll=[coll;da];

%end of the loop for search of the daughters from the previously identified mothers

%add here the daughter IDs to the tmp table line?xxx

end

%combine the collected values into the final table
	
gen_table=[gen_table;tmp];

%liberate the temporary array for the next round
tmp=[];

%redefine the newly found offspring into the search for the next round

st=reshape(coll,1,[]);

%remove cells that have no further offspring
st(isnan(st))=[]; 

gen=gen+1; %increase the generation counter

end %end of first while loop to go through all cells until no more mothers

% then add the generation times per cell in a new column (6)

for jj=1:size(gen_table,1)

	if ~isnan(gen_table(jj,5))
	gen_time=(gen_table(jj,5)-gen_table(jj,4))/time_freq; %in hours
	gen_table(jj,6)=gen_time; %in hours

	else
	gen_table(jj,6)=nan;

	end

end

Variables={'generation', 'ID_birth', 'ID_death', 'cell_birth_tp', 'cell_death_tp','generation_time'};
gen_table=array2table(gen_table,'VariableNames',Variables);

%update the mother table

sub=renamevars(sub,'MotherID','Mother_ID');

% convert generation table to time by all cell IDs,
% shape by time_length

time_length=(1:max(sub.Timepoint));

mother_ID_table=zeros(size(cellcounts.lineage(k).gen_table,1),size(time_length,2));
mother_ID_table(mother_ID_table==0)=nan;

% now fill this successively with each mother ID and its path to the next division
% for this we need the sub table with all intermediate cell IDs

lin=gen_table;

for i=1:size(lin,1)

	ID_birth=[]; ID_death=[]; dID=[]; r=[]; r2=[]; aa=[]; hlp=[];  tmp=[];
    cell_birth_tp=[]; cell_death_tp=[];
	ID_birth=lin.ID_birth(i);
	ID_death=lin.ID_death(i);

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

%remove rows with all nan

mother_ID_table(~any(~isnan(mother_ID_table), 2),:)=[];

%remove values that are zero

mother_ID_table(mother_ID_table==0)=nan;

end
end

% we have now updated the Lineages cell array, which is used for the next step

%% Section three: collect the final genealogy tables in the cellcount structure
% to make a genealogy tree

clearvars -except Lineages

%take the data belonging to a single lineage and go through all

for k=1:size(Lineages,2) %

sub=Lineages{k};

%what is the time frequency (in per hour) of the images

time_freq=3;

st=sub.MotherID(1); %st is start

gen_table=[]; %to collect the data for the generation table
gen=1; %generation counter

%go through the whole table and stop when st becomes empty
%this is the sign that all cells have been attributed

while ~isempty(st)

%to collect the daughter data from the mother cell(s) for the next round of searches

coll=[]; 

%start the loop to find new daughters from the defined mothers in this round

for j=1:length(st)

da=[]; %individual daughter(s) from each single mother
ro=[]; cell_birth=[]; cell_death=[];
ID_birth=st(j); ID_death=[];

%r is the row number of the mother cell
r=find(sub.MotherID==st(j));

%to collect the connecting rows from mother to its divided daughters, 
%which we need to define the point of cell death in the end

cell_rows=[]; 

%verify if a connecting cell exists and then connect
%if it doesn't - then break the loop and add only the cell birth values

	if ~isempty(r)

		%record cell birth as where this cell first appeared in the daughter columns
		[ro,co]=find(sub.DaughterID1==st(j) | sub.DaughterID2==st(j));

			if isempty(ro)
				%this situation arises only for the first founder cell
				cell_birth=sub.Timepoint(r);  
			else
				%sometimes happens that cell is found twice by error
				%then take the earliest timepoint
				if length(ro)>1
					cell_birth=sub.Timepoint(min(ro));
				else
					cell_birth=sub.Timepoint(ro);
				end
			end

		da=[sub.DaughterID1(r),sub.DaughterID2(r)]; 


%now there are only four possibilities: 
%either the cells divide; in that case there will be two da-values
%or there is only a single connecting cell: no division yet
%then the value of the second daughterID is NaN
%or there are no connecting cells
%then both values are NaN
%or there is no new row to be found - end of table reached

		if ~isnan(da(2)) %immediate new daughter cells
			cell_rows=[r];

			%take the same timepoint where this cell divided as its cell_death
			cell_death=sub.Timepoint(cell_rows); 

			%take the ID of the same timepoint before the division
			%ID_death=sub.DaughterID1(cell_rows);  
			ID_death=sub.MotherID(cell_rows);

			%collect the values in a temporary array
			tmp(j,:)=[gen,ID_birth,ID_death,cell_birth,cell_death];


		elseif ~isnan(da(1)) %now we have to continue to search
			while isnan(da(2))

			r=find(sub.MotherID==da(1)); %the row number of the connecting cell
			da=[sub.DaughterID1(r),sub.DaughterID2(r)]; 

			%the new value(s) of the daughter(s)
			cell_rows=[cell_rows;r];

				%if the next round finds no more mothers
				%because end of table, stop

				if isempty(r)	
					break
				end

				%sometimes there can be an error in the list
				%when the mother has the same ID as its daughter
				%this would cause an endless loop, 
				%so this has to be detected before continuing

				if da(1)==sub.MotherID(r)
					da(1)=NaN;
				elseif da(2)==sub.MotherID(r)
					da(2)=NaN;	
				end


				%if the next round has two NaN, stop

				if isempty(find(sub.MotherID==da(1), 1))
					break
				end


			end %end of the while loop

			%another check for the case where the immediate next cell will
			%have no daughters and the cell_rows would only be length 1
			%then we cannot establish the cell death
			%we add one number artifically, that can be subtracted in the
			%following line

			if isempty(cell_rows)

			%take the last timepoint of the table as its cell_death
			cell_death=max(sub.Timepoint); 

			%take the ID of the last timepoint before the division
			ID_death=sub.DaughterID1(find(sub.MotherID==st(j)));  

			%collect the data in the tmp array
			tmp(j,:)=[gen,ID_birth,ID_death,cell_birth,cell_death];

			%
			elseif length(cell_rows)==1

			cell_rows=[cell_rows;1];
	
			%take the ID of the last timepoint before the division
			ID_death=sub.MotherID(cell_rows(end-1));  

			%take the last timepoint before this cell divided as its cell_death
			cell_death=min(sub.Timepoint(find(sub.DaughterID1==ID_death))); 

			%collect the data in the tmp array
			tmp(j,:)=[gen,ID_birth,ID_death,cell_birth,cell_death];

			%
			elseif isnan(da(1)) && isnan(da(2))

				if ~isempty(cell_rows)

				cell_rows=cell_rows;
				%take the last timepoint before this cell divided as its cell_death
				cell_death=sub.Timepoint(cell_rows(end-1)); 

				%take the ID of the last timepoint before the division
				ID_death=sub.DaughterID1(cell_rows(end-1));  

				%collect the data in the tmp array
				tmp(j,:)=[gen,ID_birth,ID_death,cell_birth,cell_death];


				else
				cell_rows=r;

				%take the same timepoint where this cell divided as its cell_death
				cell_death=sub.Timepoint(cell_rows); 

				%take the ID of the same timepoint before the division
				ID_death=sub.DaughterID1(cell_rows);  

				%collect the values in a temporary array
				tmp(j,:)=[gen,ID_birth,NaN,cell_birth,NaN];

				end

			%
			else

			%take the last timepoint before this cell divided as its cell_death
			cell_death=sub.Timepoint(cell_rows(end-1)); 

			%take the ID of the last timepoint before the division
			ID_death=sub.DaughterID1(cell_rows(end-1));  

			%collect the data in the tmp array
			tmp(j,:)=[gen,ID_birth,ID_death,cell_birth,cell_death];

			end

		else %both values are NaN

			cell_rows=[r];

			%take the same timepoint where this cell divided as its cell_death
			cell_death=sub.Timepoint(cell_rows); 

			%take the ID of the same timepoint before the division
			ID_death=sub.DaughterID1(cell_rows);  

			%collect the values in a temporary array
			tmp(j,:)=[gen,ID_birth,NaN,cell_birth,NaN];

		
		end

	else 
	
	%no connecting cell has been found, we assume it died or disappeared
	%its ID_death and time of cell_death will become NaN
	% cell_birth is where it appeared

	tmp(j,:)=[gen,ID_birth,NaN,sub.Timepoint(sub.DaughterID1==ID_birth | sub.DaughterID2==ID_birth),NaN];
	
	end

	%collecting the newly found daugthers in the loop for the next round

	coll=[coll;da];

%end of the loop for search of the daughters from the previously identified mothers

%add here the daughter IDs to the tmp table line?xxx

end

%combine the collected values into the final table
	
gen_table=[gen_table;tmp];

%liberate the temporary array for the next round
tmp=[];

%redefine the newly found offspring into the search for the next round

st=reshape(coll,1,[]);

%remove cells that have no further offspring
st(isnan(st))=[]; 

gen=gen+1; %increase the generation counter

end %end of first while loop to go through all cells until no more mothers

% then add the generation times per cell in a new column (6)

for jj=1:size(gen_table,1)

	if ~isnan(gen_table(jj,5))
	gen_time=(gen_table(jj,5)-gen_table(jj,4))/time_freq; %in hours
	gen_table(jj,6)=gen_time; %in hours

	else
	gen_table(jj,6)=nan;

	end

end


%%collect in the final cellcounts structure array

Variables={'generation', 'ID_birth', 'ID_death', 'cell_birth_tp', 'cell_death_tp','generation_time'};
t=array2table(gen_table,'VariableNames',Variables);

cellcounts.lineage(k).gen_table=t;


clearvars sub gen_table gen
 
end %end of lineages loop (k)

save Lineages.mat Lineages
save cellcounts.mat cellcounts
