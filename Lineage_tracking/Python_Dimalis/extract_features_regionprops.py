import cv2
import numpy as np
import pandas as pd
from skimage import io
import sys
import glob
import dask
from helper import compute_centroids

# scipy package to compute distance between matrices of cell centroid coordinates
from scipy.spatial.distance import cdist
# skimage package to identify objects, export region properties
from skimage.measure import label, regionprops, regionprops_table

#Creation of the feature files found in feature_tables

@dask.delayed
def _compute(files_list, tp, output_dir):
    #print('Processing file ', files_list[tp])
    img1 = io.imread(files_list[tp])
    
    # See how many cells were identified in the images (ranged in the "unique" vectors)
    unique1, counts1 = np.unique(img1, return_counts=True)
    
    ###################################################
    # compute the mask centroids of all cells in img1 #
    ###################################################
    my_array1 = compute_centroids(unique1, img1)
    
    ###################################################################################
    # Compute distances between cells (in a similar way to R's pdist::pdist function) #
    ###################################################################################
    dist_mat = cdist(my_array1, my_array1)
    
    if (dist_mat.shape[0] <= 3):
        first_closest = np.repeat(0, dist_mat.shape[0])
        second_closest = np.repeat(0, dist_mat.shape[0])
    else:
        # extract 1st and 2nd closest distance of each cell's neighbours
        first_closest = [sorted(i)[1] for i in zip(*dist_mat)]
        second_closest = [sorted(i)[2] for i in zip(*dist_mat)]
    
    ###############################################################################
    # extract features for all cells using the skimage regionprops table function #
    ###############################################################################
    
    props = regionprops_table(img1, properties=['area', 'bbox', 'centroid', 'eccentricity', 'euler_number', 'extent', 'axis_minor_length', 'axis_major_length', 'feret_diameter_max', 'orientation', 'perimeter', 'solidity'])
    
    # turn props into a panda data.frame
    data_tmp = pd.DataFrame(props)
    
    # add closest neighbour distance columns
    data_tmp['Neighbors_FirstClosestDistance'] = first_closest
    data_tmp['Neighbors_SecondClosestDistance'] = second_closest
    
    #########################################################
    # Export feature cell information for the current image #
    #########################################################
    data_tmp.to_csv(output_dir + f"feature_table_time{tp}.csv")



def _extract_features(files_list, output_dir):
    delayed_results = []
    for tp in range(0,len(files_list)):
        delayed_results.append(_compute(files_list, tp, output_dir))
    dask.compute(delayed_results, scheduler='threads')

if __name__=="__main__":
    files_list = sorted(glob.glob("*.tif"))
    output_dir = sys.argv[1]
    _extract_features(files_list, output_dir)
