import pandas as pd
import os
import glob
import sys
from functools import reduce # for merging multiple dataframes
import re

#output_dir = sys.argv[1]

script_dir = os.path.dirname(os.path.abspath(__file__))
output_dir = script_dir
# load the strack cell table
strack_file = '/STrack/tracked_cells_table.xlsx'
strack_table = pd.read_excel(f"{output_dir}{strack_file}")

###################################################################
# import the feature tables and concatenate them for all timepoints
#      format them and finally merge with the strack table
###################################################################
os.chdir(f"{output_dir}{'/feature_tables/'}")

# List files in current folder
#files_list = sorted(glob.glob("*.csv"), key=len) # order files by time of creation

# List all .csv files in the current directory
files_list = glob.glob("*.csv")

# Sort the files based on the number before '_csv' in their filename
files_list = sorted(files_list, key=lambda x: int(re.search(r'_time(\d+)', x).group(1)))
# merge tables into one big panda dataframe
ft_df_list = []
#ft_df_append = pd.DataFrame()
unique_time = strack_table['Timepoint'].unique()
unique_time = sorted(unique_time)
k = 0
i = unique_time[k] #0
#append all files together
for file in files_list:
    ft_df_temp = pd.read_csv(file)
    # add a column corresponding to the timepoint
    ft_df_temp2 = ft_df_temp.assign(Timepoint=[k] * len(ft_df_temp.index))
    ft_df_list.append(ft_df_temp2)
   # ft_df_append = ft_df_append.append(ft_df_temp2, ignore_index=True)
    k += 1
    #i = unique_time[k]

ft_df_append = pd.concat(ft_df_list, ignore_index=True)
# remove columns that we don't need
ft_df_append.drop(columns=["Unnamed: 0", "bbox-0", "bbox-1", "bbox-2", "bbox-3", "euler_number"], inplace=True)

# rename columns so that they match the columns in strack table
ft_final_df = ft_df_append.rename(columns={'centroid-1': 'Centroid_x', 'centroid-0': 'Centroid_y'})

#######Attempt to duplicate columns

ft_final_df['Centroid_x_copy'] = ft_final_df['Centroid_x']  # Duplicate 'col1'
ft_final_df['Centroid_y_copy'] = ft_final_df['Centroid_y']  # Duplicate 'col2'

# truncate x and y coordinates in ft table so that they match the ones in strack table
ft_final_df.Centroid_x = [int(item) for item in ft_final_df.Centroid_x] #Here truncated but why?
ft_final_df.Centroid_y = [int(item) for item in ft_final_df.Centroid_y]

# merge strack table and ft tables (appended)
# doc found here: https://pandas.pydata.org/docs/user_guide/merging.html
result = pd.merge(strack_table, ft_final_df, how="outer", on=["Centroid_x", "Centroid_y", "Timepoint"]) #Merge on the centroids 

result['Centroid_x'] = result['Centroid_x_copy'].combine_first(result['Centroid_x']) #Replace values of centroid by the non rounded if there is no nan value
result['Centroid_y'] = result['Centroid_y_copy'].combine_first(result['Centroid_y']) #Replace values of centroid by the non rounded if there is no nan value

################################################################
# import the fluo tables and concatenate them for all timepoints
# format them and finally merge with the strack + ft merged table
################################################################

def format_fluo_tables(tmp_dir):
    # import fluo tables
    os.chdir(tmp_dir)

    # List files in current folder
    fluo_files_list = sorted(glob.glob("*.csv"), key=len) # order files by time of creation

    # merge tables into one big panda dataframe
    #fluo_df_append = pd.DataFrame()
    fluo_df_list = []
    i = 0
    #append all files together
    for file in fluo_files_list:
        fluo_df_temp = pd.read_csv(file)
        # add a column corresponding to the timepoint
        fluo_df_temp2 = fluo_df_temp.assign(Timepoint=[i] * len(fluo_df_temp.index))
        fluo_df_list.append(fluo_df_temp2)
        #fluo_df_append = fluo_df_append.append(fluo_df_temp2, ignore_index=True)
        i += 1
        
    fluo_df_append = pd.concat(fluo_df_list, ignore_index=True)
    # remove columns that we don't need and the ones that might cause confusion when concatenating with other dfs
    fluo_df_append.drop(columns=["Unnamed: 0", "Mask_nb"], inplace=True)
    # truncate x and y coordinates in fluo table so that they match the ones in strack table
    fluo_df_append.Centroid_x = [int(item) for item in fluo_df_append.Centroid_x]
    fluo_df_append.Centroid_y = [int(item) for item in fluo_df_append.Centroid_y]

    return fluo_df_append

# import and format fluo tables for the GFP channel
GFP_table = format_fluo_tables(f"{output_dir}{'/fluo_channels/GFP_results'}")
# import and format fluo tables for the mCherry channel
mCherry_table = format_fluo_tables(f"{output_dir}{'/fluo_channels/mCherry_results'}")

# merge all fluorescent tables into ones
fluo_data_frames = [GFP_table, mCherry_table]
fluo_merged = reduce(lambda  left,right: pd.merge(left,right,on=["Centroid_x", "Centroid_y", "Timepoint"],
                                            how='outer'), fluo_data_frames)

# merge result table (containing strack table + ft tables with fluo tables
result2 = pd.merge(result, fluo_merged, how="outer", on=["Centroid_x", "Centroid_y", "Timepoint"])
# remove rows where Mask_nb = NA
result3 = result2.dropna(subset=['Mask_nb'])

# Export the final panda dataframe in xlsx format
result3.to_excel(output_dir + "/final_merged_table.xlsx", index=False)
