import pandas as pd
import os
import glob
import sys

## previous script relied only on x,y coordinates to decide whether a new cell deserved a new mask_nb or not
## ISSUE: in very crowded images, a cell could appear many timepoints later with the exact same coordinates, and would get the smae mask id, meading to wrong tracks across timepoints
## solution: now, stocked_IDs has x, y coords but ALSO timepoint stored, to prevent this from happening

#Creation of tracked_cells_table.xlsx

output_dir = sys.argv[1]
os.chdir(output_dir)

# List files in current folder
files_list = sorted(glob.glob("*.csv"), key=len) # order files by numeric order
#print("attempt to debug, current path is ")
#print(os.getcwd())
#print(files_list)
# merge tables into one big panda dataframe
df_append = pd.DataFrame()
#append all files together
for file in files_list:
    df_temp = pd.read_csv(file)
    df_append = df_append.append(df_temp, ignore_index=True)

# loop over the tables rows and modify mother and daughter mask IDs
# such that they are fit for building a graph
stocked_IDs = {}
incremental_ID = 1

for i, row in df_append.iterrows():
    # change mother iD
    coords_tmp = f"{row.Centroid_x_mother},{row.Centroid_y_mother},{row.Timepoint-1}"
    if coords_tmp not in stocked_IDs.keys():
        stocked_IDs[coords_tmp] = incremental_ID
        incremental_ID += 1
    df_append['Mother_mask'][i] = stocked_IDs[coords_tmp]

    # change daughter ID
    coords_tmp_d = f"{row.Centroid_x},{row.Centroid_y},{row.Timepoint}"
    if coords_tmp_d not in stocked_IDs.keys():
        stocked_IDs[coords_tmp_d] = incremental_ID
        incremental_ID += 1
    df_append['Mask_nb'][i] = stocked_IDs[coords_tmp_d]

# Finally, remove the redundant rows, where Mother_mask == Mask_nb
same_value_rows = df_append[df_append['Mother_mask'] == df_append['Mask_nb']]
df_append = df_append.drop(same_value_rows.index)
# and reset index: drop=True ensures that the old index is dropped and a new default index is assigned to the DataFrame.
df_append = df_append.reset_index(drop=True)

# Export the final panda dataframe in csv format
df_append.to_excel(output_dir + "complete_tracking_table.xlsx", index=False)

# Select the mother and daughter cells info
daughter_cells = df_append[["Mask_nb", "Centroid_x", "Centroid_y", "Timepoint"]]
mother_cells = df_append[["Mother_mask", "Centroid_x_mother", "Centroid_y_mother", "Timepoint"]]
# remove one to the mother cell timepoints
mother_cells.loc[:, "Timepoint"] = mother_cells["Timepoint"].apply(lambda x: x - 1)

# rename the columns in the mother cells dataframe to match the daughter_cells df
mother_cells = mother_cells.rename(columns={"Mother_mask": "Mask_nb", "Centroid_x_mother": "Centroid_x", "Centroid_y_mother": "Centroid_y"})

# Merge both tables, and remove duplicates:
all_cells = pd.concat([mother_cells,daughter_cells]).drop_duplicates().reset_index(drop=True)

# Export the panda dataframe containing info on all cells in excel format
all_cells.to_excel(output_dir + "tracked_cells_table.xlsx", index=False)
