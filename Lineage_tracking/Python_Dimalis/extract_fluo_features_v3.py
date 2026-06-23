import cv2
import numpy as np
import pandas as pd
from skimage import io
import os # to list subfolders in a specified directory, create folders ...
import sys
import glob
# skimage package to identify objects, export region properties
from skimage.measure import label, regionprops, regionprops_table
import dask


#Creation fluo channel files .csv, found in results

@dask.delayed
def _compute(tp, segmented_masks_list, tmp_files_list, input_dir):
    # Add a 0 in front of tp if contains only one number
    tp1 = '%02d' % (tp)
    # Import segmented masks corresponding to timepoint tp
    #print('Processing file for timepoint ', tp1)
    img0 = io.imread(segmented_masks_list[tp])
    # Import corresponding fluorescence images
    img0_fluo = io.imread(tmp_files_list[tp])

    # count cell masks in the current image
    unique0, counts0 = np.unique(img0, return_counts=True)

    # Create empty result table
    fluo_results = pd.DataFrame(columns=['Timepoint', 'Mask_nb', 'Centroid_x', 'Centroid_y', fluo_tmp_name + '_mean', fluo_tmp_name + '_median', fluo_tmp_name + '_1st_quartile', fluo_tmp_name +  '_3rd_quartile', fluo_tmp_name + '_sum', fluo_tmp_name + '_min', fluo_tmp_name + '_max'])
    for mask_tmp in unique0[1:]:
        # Keep only current mask, replace all other values by 0
        img_tmp = np.where(img0 != mask_tmp, 0, img0)
        # Keep fluorescence only in the current mask, replace all other fluorescence by 0
        img_tmp_fluo = np.where(img0 != mask_tmp, 0, img0_fluo)

        # Measure total, average, 1st quartile, 3rd quartile, min and max fluorescence in the current cell mask
        fluo_mean = img_tmp_fluo[np.nonzero(img_tmp_fluo)].mean()
        fluo_median = np.quantile(img_tmp_fluo[np.nonzero(img_tmp_fluo)], 0.5)
        fluo_1quart = np.quantile(img_tmp_fluo[np.nonzero(img_tmp_fluo)], 0.25)
        fluo_3quart = np.quantile(img_tmp_fluo[np.nonzero(img_tmp_fluo)], 0.75)
        fluo_sum = img_tmp_fluo[np.nonzero(img_tmp_fluo)].sum()
        fluo_min = img_tmp_fluo[np.nonzero(img_tmp_fluo)].min()
        fluo_max = img_tmp_fluo[np.nonzero(img_tmp_fluo)].max()

        # compute cell centroid
        M = cv2.moments(img_tmp)
        cx = int(M['m10'] / M['m00'])
        cy = int(M['m01'] / M['m00'])

        # save results in a panda dataframe
        fluo_results.loc[mask_tmp] = [tp, # timepoint
                                    mask_tmp, # mask number
                                    cx, # = cx
                                    cy, # = cy
                                    fluo_mean, fluo_median, fluo_1quart, fluo_3quart,
                                    fluo_sum, fluo_min, fluo_max]

        # Save the panda dataframe for the current image timepoint
        fluo_results.to_csv(input_dir + "/fluo_channels/" + fluo_tmp_name + "_results/" + fluo_tmp_name + f"_table_time{tp}.csv")


def _extract_fluo_features(segmented_masks_list, tmp_files_list, input_dir):
    delayed_results = []
    for tp in range(len(segmented_masks_list)):
        delayed_results.append(_compute(tp, segmented_masks_list, tmp_files_list, input_dir))
    dask.compute(delayed_results, scheduler='threads')

if __name__=="__main__":
    input_dir = sys.argv[1]
    # list files in Omnipose folder
    segmented_masks_list = sorted(glob.glob(input_dir + "/Omnipose/29.5/"+ "*.tif"))
    # !! mac can create .DS_store folders that break the code: I select only folders that end with "channel"
    included_extensions = ['channel']
    fluo_channels_list = [fn for fn in os.listdir(input_dir + '/fluo_channels/')
    if any(fn.endswith(ext) for ext in included_extensions)]
    print("fluo channels are ", fluo_channels_list)

    for fluo_tmp in fluo_channels_list:
        #print("Analysing fluo results for ", fluo_tmp)
        fluo_tmp_name = fluo_tmp.replace("_channel", "")
    
        # create subfolder to output the result tables into
        os.makedirs(input_dir + "/fluo_channels/" + fluo_tmp_name + "_results/", exist_ok=True)
    
        # list files in the current fluorescence folder
        folder_to_list_files_in = input_dir + "/fluo_channels/" + fluo_tmp + "/"
        tmp_files_list = sorted(glob.glob(folder_to_list_files_in + "*.tif"))
        _extract_fluo_features(segmented_masks_list, tmp_files_list, input_dir)
