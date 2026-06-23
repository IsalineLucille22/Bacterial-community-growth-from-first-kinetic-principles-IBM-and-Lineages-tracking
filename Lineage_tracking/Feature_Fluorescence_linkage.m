%% This script links the features to the various fluorescence data from Bouke's
% experiments.
% go to the main position folder to start the linkage and conversion
% an output is saved in form of a 'final_merged_table.xlsx' that can
% be used in the cell linkage script

feature_fds = fileDatastore('feature_tables/feature_table*.csv', 'ReadFcn', @importdata);
feature_GFP = fileDatastore('fluo_channels/GFP_results/GFP_table*.csv', 'ReadFcn', @importdata);
feature_mCherry = fileDatastore('fluo_channels/mCherry_results/mCherry_table*.csv', 'ReadFcn', @importdata);
feature_CFP = fileDatastore('fluo_channels/CFP_results/CFP_table*.csv', 'ReadFcn', @importdata);
feature_YFP = fileDatastore('fluo_channels/YFP_results/YFP_table*.csv', 'ReadFcn', @importdata);
feature_DAPI = fileDatastore('fluo_channels/DAPI_results/DAPI_table*.csv', 'ReadFcn', @importdata);


fullFileNames_fds = feature_fds.Files;
fullFileNames_GFP = feature_GFP.Files;
fullFileNames_mCherry = feature_mCherry.Files;
fullFileNames_CFP = feature_CFP.Files;
fullFileNames_YFP = feature_YFP.Files;
fullFileNames_DAPI = feature_DAPI.Files;

tracked = readtable('Strack/tracked_cells_table.xlsx'); %to have the correct Mask_nb
times=unique(tracked.Timepoint);

numFiles=length(fullFileNames_fds);

total_features=[]; %to collect the data

for i=1:numFiles

feature_table = readtable(string(fullFileNames_fds(i)));
feature_table = removevars(feature_table,{'Var1'});

GFP_table = readtable(string(fullFileNames_GFP(i)));
GFP_table = removevars(GFP_table,{'Var1'});

X=[GFP_table.Centroid_x,GFP_table.Centroid_y];
Y=[feature_table.centroid_1,feature_table.centroid_0];
[E,I]=pdist2(X,Y,'Euclidean','smallest',1);

GFP_table=GFP_table(I,:); %sort according to the distance

mCherry_table=readtable(string(fullFileNames_mCherry(i)));
mCherry_table = removevars(mCherry_table,{'Var1'});

X=[mCherry_table.Centroid_x,mCherry_table.Centroid_y];
Y=[feature_table.centroid_1,feature_table.centroid_0];
[E,I]=pdist2(X,Y,'Euclidean','smallest',1);

mCherry_table=mCherry_table(I,:);
mCherry_table = removevars(mCherry_table,{'Mask_nb','Centroid_x','Centroid_y','Timepoint'});

YFP_table=readtable(string(fullFileNames_YFP(i)));
YFP_table = removevars(YFP_table,{'Var1'});

X=[YFP_table.Centroid_x,YFP_table.Centroid_y];
Y=[feature_table.centroid_1,feature_table.centroid_0];
[E,I]=pdist2(X,Y,'Euclidean','smallest',1);

YFP_table=YFP_table(I,:);
YFP_table = removevars(YFP_table,{'Mask_nb','Centroid_x','Centroid_y','Timepoint'});

CFP_table=readtable(string(fullFileNames_CFP(i)));
CFP_table = removevars(CFP_table,{'Var1'});

X=[CFP_table.Centroid_x,CFP_table.Centroid_y];
Y=[feature_table.centroid_1,feature_table.centroid_0];
[E,I]=pdist2(X,Y,'Euclidean','smallest',1);

CFP_table=CFP_table(I,:);
CFP_table = removevars(CFP_table,{'Mask_nb','Centroid_x','Centroid_y','Timepoint'});

DAPI_table=readtable(string(fullFileNames_DAPI(i)));
DAPI_table = removevars(DAPI_table,{'Var1'});

X=[DAPI_table.Centroid_x,DAPI_table.Centroid_y];
Y=[feature_table.centroid_1,feature_table.centroid_0];
[E,I]=pdist2(X,Y,'Euclidean','smallest',1);

DAPI_table=DAPI_table(I,:);
DAPI_table = removevars(DAPI_table,{'Mask_nb','Centroid_x','Centroid_y','Timepoint'});

time_track=tracked(tracked.Timepoint==times(i),:);
X=[time_track.Centroid_x, time_track.Centroid_y];
Y=[feature_table.centroid_1,feature_table.centroid_0];
[E,I]=pdist2(X,Y,'Euclidean','smallest',1);

time_track=time_track(I,:);
GFP_table=removevars(GFP_table,{'Timepoint','Mask_nb','Centroid_x','Centroid_y'});

fused_table=[time_track,feature_table,GFP_table,mCherry_table,CFP_table,YFP_table,DAPI_table];

total_features=vertcat(total_features,fused_table);

clearvars -except total_features fullFileNames_fds i fullFileNames_GFP fullFileNames_CFP fullFileNames_DAPI fullFileNames_mCherry fullFileNames_YFP tracked times

end


writetable(total_features, 'final_merged_table.xlsx');